"""Reader for Agilent MassHunter GC/MS ``.d`` datasets (MSScan.bin + MSPeak.bin).

Pure standard library -- no numpy/pandas required.

Binary layout
-------------
``AcqData/MSScan.bin`` is an int32 at offset 0x58 giving the start of a packed
array of fixed-size scan records.  The record fields are exactly the sequence
declared in ``AcqData/MSScan.xsd`` (ScanRecordType followed by one
SpectrumParamsType), stored little-endian with no padding -> 184 bytes.

``AcqData/MSPeak.bin`` holds the centroid peak lists.  For each scan,
``SpectrumOffset`` points at ``PointCount`` little-endian float32 m/z values
immediately followed by ``PointCount`` float32 abundances (``ByteCount`` ==
``PointCount`` * 8).

This layout was validated against every scan of this dataset: the sum of the
decoded abundances reproduces the stored ``TIC``, the abundance argmax
reproduces ``BasePeakMZ``/``BasePeakValue``, the m/z arrays are monotonic and
span ``MinX``..``MaxX``, and the last peak block ends exactly at EOF.

Only the standard library is required.  If numpy happens to be installed,
:meth:`AgilentDotD.flat_peaks` and :meth:`AgilentDotD.eic` become available and
give fast vectorised extracted-ion chromatograms; :meth:`AgilentDotD.xic` is the
pure-python equivalent.
"""

from __future__ import annotations

import base64
import hashlib
import mmap
import os
import struct
import xml.etree.ElementTree as ET
import zlib
from typing import Dict, Iterator, List, NamedTuple, Optional, Sequence, Tuple

__all__ = ["AgilentDotD", "ScanRecord", "find_datasets", "EIC_EDGE"]

# ---------------------------------------------------------------------------
# MSScan.bin record
# ---------------------------------------------------------------------------

_REC = struct.Struct("<3i d 2i 3d 4i 7d 2i i q 2i 4d")
assert _REC.size == 184
_DATA_OFFSET_POS = 0x58

#: Window slack (u) added on each side of an EIC target -- see
#: :meth:`AgilentDotD.eic_window`.  Far below the 0.1 u acquisition step.
EIC_EDGE = 1e-4


class ScanRecord(NamedTuple):
    scan_id: int
    scan_method_id: int
    time_segment_id: int
    scan_time: float           # minutes
    ms_level: int
    scan_type: int
    tic: float
    base_peak_mz: float
    base_peak_value: float
    cycle_number: int
    status: int
    ion_mode: int
    ion_polarity: int
    fragmentor: float
    collision_energy: float
    mz_of_interest: float
    sampling_period: float
    mass_range_min: float
    mass_range_max: float
    threshold: float
    is_fragmentor_dynamic: int
    is_collision_energy_dynamic: int
    spectrum_format_id: int
    spectrum_offset: int
    byte_count: int
    point_count: int
    min_y: float
    max_y: float
    min_x: float
    max_x: float


# Agilent enum -> label.  Values seen in this dataset are marked; the rest come
# from the MassHunter data-access documentation and are best-effort labels.
SCAN_TYPES = {1: "Scan", 2: "SelectedIon", 4: "HighResolutionScan",
              8: "TotalIon", 256: "MultipleReaction", 512: "ProductIon",
              1024: "PrecursorIon", 2048: "NeutralGain", 4096: "NeutralLoss"}
ION_MODES = {0: "Unspecified", 1: "Mixed", 2: "EI", 4: "CI", 8: "Maldi",
             16: "Appi", 32: "Apci", 64: "Esi", 128: "NanoEsi",
             512: "MsChip", 1024: "ICP", 2048: "Jetstream"}
ION_POLARITY = {0: "positive", 1: "negative", 2: "unassigned", 3: "mixed"}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _xml_root(path: str):
    if not os.path.exists(path):
        return None
    # Agilent writes some of these with a UTF-8 BOM.
    with open(path, "rb") as fh:
        raw = fh.read()
    if raw.startswith(b"\xef\xbb\xbf"):
        raw = raw[3:]
    try:
        return ET.fromstring(raw)
    except ET.ParseError:
        return None


def _localname(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def decode_text_file(path: str) -> Optional[str]:
    """Decode an Agilent text report (UTF-16LE with BOM, or UTF-8)."""
    if not os.path.exists(path):
        return None
    with open(path, "rb") as fh:
        raw = fh.read()
    for bom, enc in ((b"\xff\xfe", "utf-16-le"), (b"\xfe\xff", "utf-16-be"),
                     (b"\xef\xbb\xbf", "utf-8")):
        if raw.startswith(bom):
            return raw[len(bom):].decode(enc, "replace")
    for enc in ("utf-8", "gb18030", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", "replace")


def find_datasets(root: str) -> List[str]:
    """All ``*.d`` directories under *root* (recursively), sorted."""
    out = []
    for dirpath, dirnames, _ in os.walk(root):
        for d in list(dirnames):
            if d.lower().endswith(".d"):
                out.append(os.path.join(dirpath, d))
                dirnames.remove(d)          # don't descend into the .d
    return sorted(out)


# ---------------------------------------------------------------------------
# main reader
# ---------------------------------------------------------------------------

class AgilentDotD:
    def __init__(self, path: str):
        self.path = os.path.abspath(path.rstrip(os.sep))
        self.name = os.path.basename(self.path)
        if self.name.lower().endswith(".d"):
            self.name = self.name[:-2]
        self.acqdata = os.path.join(self.path, "AcqData")
        if not os.path.isdir(self.acqdata):
            raise ValueError("not an Agilent .d dataset (no AcqData): %s" % path)

        self._scans: Optional[List[ScanRecord]] = None
        self._flat = None
        self._peak_mm = None
        self._peak_fh = None

        self.sample_info = self._read_sample_info()
        self.contents = self._read_contents()
        self.devices = self._read_devices()
        self.time_segments = self._read_time_segments()
        self.method_dir = self._find_method_dir()
        self.ms_method = self._read_ms_method()

    # -- context manager ---------------------------------------------------
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def close(self):
        if self._peak_mm is not None:
            self._peak_mm.close()
            self._peak_mm = None
        if self._peak_fh is not None:
            self._peak_fh.close()
            self._peak_fh = None

    # -- metadata ----------------------------------------------------------
    def _read_sample_info(self) -> Dict[str, str]:
        root = _xml_root(os.path.join(self.acqdata, "sample_info.xml"))
        info: Dict[str, str] = {}
        if root is None:
            return info
        for field in root.findall("Field"):
            name = (field.findtext("Name") or "").strip()
            disp = (field.findtext("DisplayName") or "").strip()
            val = (field.findtext("Value") or "").strip()
            # Some rows have a numeric hash as <Name>; the DisplayName is the
            # only human-readable key in that case.
            key = name if name and not name.isdigit() else disp
            if key:
                info[key] = val
        return info

    def _read_contents(self) -> Dict[str, str]:
        root = _xml_root(os.path.join(self.acqdata, "Contents.xml"))
        if root is None:
            return {}
        return {_localname(c.tag): (c.text or "").strip() for c in root}

    def _read_devices(self) -> List[Dict[str, str]]:
        root = _xml_root(os.path.join(self.acqdata, "Devices.xml"))
        if root is None:
            return []
        out = []
        for dev in root.findall("Device"):
            d = {_localname(c.tag): (c.text or "").strip() for c in dev}
            d["DeviceID"] = dev.get("DeviceID", "")
            out.append(d)
        return out

    def _read_time_segments(self) -> List[Dict[str, str]]:
        root = _xml_root(os.path.join(self.acqdata, "MSTS.xml"))
        if root is None:
            return []
        out = []
        for seg in root.findall("TimeSegment"):
            d = {_localname(c.tag): (c.text or "").strip() for c in seg}
            d["TimeSegmentID"] = seg.get("TimeSegmentID", "")
            out.append(d)
        return out

    def _find_method_dir(self) -> Optional[str]:
        for entry in sorted(os.listdir(self.acqdata)):
            p = os.path.join(self.acqdata, entry)
            if entry.lower().endswith(".m") and os.path.isdir(p):
                return p
        return None

    @property
    def method_name(self) -> str:
        if self.method_dir:
            return os.path.basename(self.method_dir)
        return self.sample_info.get("Method", "")

    def _read_ms_method(self) -> Dict[str, object]:
        """Flatten the interesting parts of ``qqqacqmethod.xml``."""
        if not self.method_dir:
            return {}
        root = _xml_root(os.path.join(self.method_dir, "qqqacqmethod.xml"))
        if root is None:
            return {}
        m: Dict[str, object] = {}
        for tag in ("msInstrument", "ionSource", "tuneFile", "stopMode",
                    "stopTime", "solventDelay", "collisionGasOn", "useGain",
                    "enableNR", "isTimeFilterEnabled", "timeFilterPeakWidth"):
            txt = root.findtext(tag)
            if txt is not None:
                m[tag] = txt.strip()
        srcs = {}
        for sp in root.iter("sourceParameter"):
            srcs[(sp.findtext("id") or "").strip()] = (
                sp.findtext("posPolarityValue") or "").strip()
        if srcs:
            m["sourceParameters"] = srcs
        segs = []
        for seg in root.iter("timeSegment"):
            s = {"index": seg.findtext("index"),
                 "startTime": seg.findtext("startTime"),
                 "scanSegments": []}
            for ss in seg.iter("scanSegment"):
                d = {"index": ss.findtext("index"),
                     "scanType": ss.findtext("scanType"),
                     "scanTime": ss.findtext("scanTime"),
                     "dataStorage": ss.findtext("dataStorage"),
                     "threshold": ss.findtext("threshold"),
                     "scanElements": []}
                for se in ss.iter("scanElement"):
                    d["scanElements"].append(
                        {_localname(c.tag): (c.text or "").strip() for c in se})
                s["scanSegments"].append(d)
            segs.append(s)
        if segs:
            m["timeSegments"] = segs
        return m

    def method_text(self) -> Optional[str]:
        """Decoded ``acqmeth.txt`` (the human-readable method printout)."""
        for cand in (os.path.join(self.path, "acqmeth.txt"),
                     os.path.join(self.method_dir or "", "acqmeth.txt")):
            if cand and os.path.exists(cand):
                return decode_text_file(cand)
        return None

    # -- convenience accessors --------------------------------------------
    @property
    def sample_name(self) -> str:
        return self.sample_info.get("Sample Name") or self.name

    @property
    def acquired_time(self) -> str:
        return (self.contents.get("AcquiredTime")
                or self.sample_info.get("AcqTime", ""))

    @property
    def instrument(self) -> str:
        return (self.contents.get("InstrumentName")
                or self.sample_info.get("InstrumentName", ""))

    @property
    def ms_device(self) -> Dict[str, str]:
        for d in self.devices:
            if d.get("Name", "").lower().startswith(("tandemquad", "quad", "ms")):
                return d
        return self.devices[0] if self.devices else {}

    # -- scan table --------------------------------------------------------
    @property
    def scans(self) -> List[ScanRecord]:
        if self._scans is None:
            self._scans = self._read_scan_table()
        return self._scans

    def _read_scan_table(self) -> List[ScanRecord]:
        p = os.path.join(self.acqdata, "MSScan.bin")
        with open(p, "rb") as fh:
            buf = fh.read()
        size = len(buf)
        if size < _DATA_OFFSET_POS + 4:
            raise ValueError("MSScan.bin too small: %s" % p)

        start = struct.unpack_from("<i", buf, _DATA_OFFSET_POS)[0]
        if not (0 < start <= size and (size - start) % _REC.size == 0):
            start = self._probe_data_offset(buf)
        n = (size - start) // _REC.size
        return [ScanRecord(*_REC.unpack_from(buf, start + i * _REC.size))
                for i in range(n)]

    @staticmethod
    def _probe_data_offset(buf: bytes) -> int:
        """Fallback: find the record array start by consistency checks."""
        size = len(buf)
        for off in range(4, min(4096, size), 4):
            if (size - off) % _REC.size:
                continue
            r = ScanRecord(*_REC.unpack_from(buf, off))
            if r.scan_id == 1 and 0.0 <= r.scan_time < 10000.0 and \
                    r.point_count >= 0 and r.byte_count == r.point_count * 8:
                return off
        raise ValueError("could not locate MSScan.bin record array")

    def __len__(self) -> int:
        return len(self.scans)

    # -- spectra -----------------------------------------------------------
    def _peaks(self) -> mmap.mmap:
        if self._peak_mm is None:
            self._peak_fh = open(os.path.join(self.acqdata, "MSPeak.bin"), "rb")
            self._peak_mm = mmap.mmap(self._peak_fh.fileno(), 0,
                                      access=mmap.ACCESS_READ)
        return self._peak_mm

    def spectrum(self, index: int) -> Tuple[Sequence[float], Sequence[float]]:
        """Centroid spectrum for scan *index* as ``(mz, abundance)``."""
        r = self.scans[index]
        if r.point_count <= 0:
            return (), ()
        mm = self._peaks()
        n = r.point_count
        mz = struct.unpack_from("<%df" % n, mm, r.spectrum_offset)
        ab = struct.unpack_from("<%df" % n, mm, r.spectrum_offset + n * 4)
        return mz, ab

    def iter_spectra(self) -> Iterator[Tuple[ScanRecord, Sequence[float],
                                             Sequence[float]]]:
        for i, r in enumerate(self.scans):
            mz, ab = self.spectrum(i)
            yield r, mz, ab

    # -- derived traces ----------------------------------------------------
    def tic(self) -> Tuple[List[float], List[float]]:
        rt = [r.scan_time for r in self.scans]
        return rt, [r.tic for r in self.scans]

    def xic(self, mz_low: float, mz_high: float) -> Tuple[List[float],
                                                          List[float]]:
        """Extracted ion chromatogram, summing peaks in [mz_low, mz_high].

        Pure python; see :meth:`eic` for the vectorised version.
        """
        rt, sig = [], []
        for r, mz, ab in self.iter_spectra():
            rt.append(r.scan_time)
            total = 0.0
            for m, a in zip(mz, ab):
                if mz_low <= m <= mz_high:
                    total += a
            sig.append(total)
        return rt, sig

    # -- numpy-accelerated extraction --------------------------------------
    def flat_peaks(self):
        """Every centroid of the run as ``(mz, abundance, scan_index)``.

        Three parallel numpy arrays covering all scans, cached on the instance.
        ``scan_index`` is the 0-based position in :attr:`scans`.  m/z and
        abundance are widened to float64 (lossless from the stored float32) so
        that window comparisons behave identically to :meth:`xic`, which works
        in python floats.
        """
        if getattr(self, "_flat", None) is not None:
            return self._flat
        try:
            import numpy as np
        except ImportError as exc:                      # pragma: no cover
            raise ImportError("flat_peaks()/eic() need numpy; use xic() "
                              "for the pure-python version") from exc

        mm = self._peaks()
        total = sum(r.point_count for r in self.scans)
        mz = np.empty(total, dtype=np.float64)
        ab = np.empty(total, dtype=np.float64)
        idx = np.empty(total, dtype=np.int32)
        pos = 0
        for i, r in enumerate(self.scans):
            n = r.point_count
            if n <= 0:
                continue
            off = r.spectrum_offset
            mz[pos:pos + n] = np.frombuffer(mm, dtype="<f4", count=n,
                                            offset=off)
            ab[pos:pos + n] = np.frombuffer(mm, dtype="<f4", count=n,
                                            offset=off + n * 4)
            idx[pos:pos + n] = i
            pos += n
        assert pos == total, (pos, total)
        self._flat = (mz, ab, idx)
        return self._flat

    def rt_array(self):
        import numpy as np
        return np.asarray([r.scan_time for r in self.scans], dtype=float)

    def eic_window(self, target: float, tol: float) -> Tuple[float, float]:
        """Inclusive m/z window used by :meth:`eic` for *target* +/- *tol*.

        Widened by :data:`EIC_EDGE` so that centroids sitting exactly on the
        window edge are included deterministically.  This matters here: the
        acquisition steps m/z in 0.1 u, values are stored as float32, and a
        request like 81 +/- 0.3 lands its edges precisely on the 80.7 and 81.3
        grid points, where float rounding would otherwise decide the outcome.
        """
        return (target - tol - EIC_EDGE, target + tol + EIC_EDGE)

    def eic(self, targets, tol: float = 0.3):
        """Extracted ion chromatograms for one or more target m/z.

        *targets* is a single m/z or an iterable of them; every centroid inside
        :meth:`eic_window` is summed into that target's trace.

        Returns ``(rt, traces)`` where ``rt`` has one entry per scan and
        ``traces`` is a 2-D array of shape ``(len(targets), n_scans)``.  A
        scalar *targets* yields a 1-D trace.
        """
        import numpy as np

        single = isinstance(targets, (int, float))
        tlist = [float(targets)] if single else [float(t) for t in targets]
        mz, ab, idx = self.flat_peaks()
        n = len(self.scans)
        out = np.zeros((len(tlist), n), dtype=float)
        for k, t in enumerate(tlist):
            lo, hi = self.eic_window(t, tol)
            sel = (mz >= lo) & (mz <= hi)
            if sel.any():
                out[k] = np.bincount(idx[sel], weights=ab[sel], minlength=n)
        rt = self.rt_array()
        return (rt, out[0]) if single else (rt, out)

    def average_spectrum(self, decimals: int = 1,
                         rt_range: Optional[Tuple[float, float]] = None
                         ) -> List[Tuple[float, float, int]]:
        """Summed spectrum binned to *decimals* places.

        Returns ``(mz, summed_abundance, n_scans_contributing)`` sorted by m/z.
        """
        acc: Dict[float, List[float]] = {}
        for r, mz, ab in self.iter_spectra():
            if rt_range and not (rt_range[0] <= r.scan_time <= rt_range[1]):
                continue
            for m, a in zip(mz, ab):
                key = round(m, decimals)
                slot = acc.get(key)
                if slot is None:
                    acc[key] = [a, 1]
                else:
                    slot[0] += a
                    slot[1] += 1
        return [(k, v[0], int(v[1])) for k, v in sorted(acc.items())]

    # -- mzML export -------------------------------------------------------
    def write_mzml(self, out_path: str, compress: bool = True) -> str:
        """Write an indexedmzML 1.1.0 file with all centroid scans + the TIC."""
        writer = _MzMLWriter(self, out_path, compress=compress)
        return writer.run()


# ---------------------------------------------------------------------------
# mzML 1.1.0 writer
# ---------------------------------------------------------------------------

def _esc(s: object) -> str:
    return (str(s).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _num(x: float) -> str:
    """Shortest decimal string that round-trips the float exactly."""
    return repr(float(x))


def _ncname(s: str, prefix: str = "run_") -> str:
    """Make *s* usable as an mzML id (XML NCName: no leading digit, no colon)."""
    out = "".join(c if (c.isalnum() or c in "_-.") else "_" for c in str(s))
    if not out or not (out[0].isalpha() or out[0] == "_"):
        out = prefix + out
    return out


def _encode_array(values: Sequence[float], compress: bool) -> Tuple[bytes, bool]:
    raw = struct.pack("<%df" % len(values), *values)
    if compress and raw:
        comp = zlib.compress(raw, 6)
        if len(comp) < len(raw):
            return base64.b64encode(comp), True
    return base64.b64encode(raw), False


class _MzMLWriter:
    def __init__(self, ds: AgilentDotD, out_path: str, compress: bool = True):
        self.ds = ds
        self.out_path = out_path
        self.compress = compress
        self.fh = None
        self.spec_offsets: List[Tuple[str, int]] = []
        self.chrom_offsets: List[Tuple[str, int]] = []

    # -- low level ---------------------------------------------------------
    def _w(self, text: str):
        self.fh.write(text.encode("utf-8"))

    def _pos(self) -> int:
        self.fh.flush()
        return self.fh.tell()

    def _binary_array(self, values: Sequence[float], kind: str, indent: str):
        """kind: 'mz' | 'intensity' | 'time'"""
        b64, compressed = _encode_array(values, self.compress)
        self._w('%s<binaryDataArray encodedLength="%d">\n' % (indent, len(b64)))
        self._w('%s  <cvParam cvRef="MS" accession="MS:1000521" '
                'name="32-bit float" value=""/>\n' % indent)
        self._w('%s  <cvParam cvRef="MS" accession="%s" name="%s" value=""/>\n'
                % (indent, "MS:1000574" if compressed else "MS:1000576",
                   "zlib compression" if compressed else "no compression"))
        if kind == "mz":
            self._w('%s  <cvParam cvRef="MS" accession="MS:1000514" '
                    'name="m/z array" value="" unitCvRef="MS" '
                    'unitAccession="MS:1000040" unitName="m/z"/>\n' % indent)
        elif kind == "intensity":
            self._w('%s  <cvParam cvRef="MS" accession="MS:1000515" '
                    'name="intensity array" value="" unitCvRef="MS" '
                    'unitAccession="MS:1000131" '
                    'unitName="number of detector counts"/>\n' % indent)
        else:
            self._w('%s  <cvParam cvRef="MS" accession="MS:1000595" '
                    'name="time array" value="" unitCvRef="UO" '
                    'unitAccession="UO:0000031" unitName="minute"/>\n' % indent)
        if b64:
            self._w('%s  <binary>%s</binary>\n' % (indent, b64.decode("ascii")))
        else:
            self._w('%s  <binary/>\n' % indent)
        self._w('%s</binaryDataArray>\n' % indent)

    # -- document ----------------------------------------------------------
    def run(self) -> str:
        ds = self.ds
        os.makedirs(os.path.dirname(os.path.abspath(self.out_path)) or ".",
                    exist_ok=True)
        with open(self.out_path, "wb") as fh:
            self.fh = fh
            self._header()
            self._spectrum_list()
            self._chromatogram_list()
            self._w("  </run>\n</mzML>\n")
            index_offset = self._pos()
            self._index()
            self._w("<indexListOffset>%d</indexListOffset>\n" % index_offset)
            self._w("<fileChecksum>")
            fh.flush()
            digest = self._sha1_upto(fh.tell())
            self._w("%s</fileChecksum>\n</indexedmzML>\n" % digest)
        self.fh = None
        return self.out_path

    def _sha1_upto(self, nbytes: int) -> str:
        h = hashlib.sha1()
        with open(self.out_path, "rb") as f:
            left = nbytes
            while left > 0:
                chunk = f.read(min(1 << 20, left))
                if not chunk:
                    break
                left -= len(chunk)
                h.update(chunk)
        return h.hexdigest()

    def _header(self):
        ds = self.ds
        dev = ds.ms_device
        model = dev.get("ModelNumber", "")
        serial = dev.get("SerialNumber", "")
        run_id = _esc(_ncname(ds.name))
        sample_name = ds.sample_name
        src_dir = os.path.basename(ds.path)

        self._w('<?xml version="1.0" encoding="utf-8"?>\n')
        self._w('<indexedmzML xmlns="http://psi.hupo.org/ms/mzml" '
                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                'xsi:schemaLocation="http://psi.hupo.org/ms/mzml '
                'http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">\n')
        self._w('<mzML xmlns="http://psi.hupo.org/ms/mzml" '
                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                'xsi:schemaLocation="http://psi.hupo.org/ms/mzml '
                'http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" '
                'version="1.1.0" id="%s">\n' % run_id)
        self._w('  <cvList count="2">\n'
                '    <cv id="MS" fullName="PSI-MS controlled vocabulary" '
                'version="4.1.0" URI="https://raw.githubusercontent.com/'
                'HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>\n'
                '    <cv id="UO" fullName="Unit Ontology" version="releases/'
                '2020-03-10" URI="http://purl.obolibrary.org/obo/uo.obo"/>\n'
                '  </cvList>\n')

        self._w('  <fileDescription>\n    <fileContent>\n'
                '      <cvParam cvRef="MS" accession="MS:1000579" '
                'name="MS1 spectrum" value=""/>\n'
                '      <cvParam cvRef="MS" accession="MS:1000127" '
                'name="centroid spectrum" value=""/>\n'
                '    </fileContent>\n')
        self._w('    <sourceFileList count="1">\n'
                '      <sourceFile id="RAW" name="%s" location="file://%s">\n'
                '        <cvParam cvRef="MS" accession="MS:1001509" '
                'name="Agilent MassHunter format" value=""/>\n'
                '        <cvParam cvRef="MS" accession="MS:1001508" '
                'name="Agilent MassHunter nativeID format" value=""/>\n'
                '      </sourceFile>\n    </sourceFileList>\n'
                % (_esc(src_dir), _esc(os.path.dirname(self.ds.path))))
        self._w('  </fileDescription>\n')

        self._w('  <sampleList count="1">\n'
                '    <sample id="SAMPLE1" name="%s">\n'
                '      <userParam name="vial" value="%s"/>\n'
                '      <userParam name="injection volume (uL)" value="%s"/>\n'
                '      <userParam name="dilution" value="%s"/>\n'
                '    </sample>\n  </sampleList>\n'
                % (_esc(sample_name),
                   _esc(ds.sample_info.get("Sample Position", "")),
                   _esc(ds.sample_info.get("Inj Vol (\u00b5l)", "")),
                   _esc(ds.sample_info.get("Dilution", ""))))

        self._w('  <softwareList count="2">\n'
                '    <software id="acquisition" version="%s">\n'
                '      <cvParam cvRef="MS" accession="MS:1000678" '
                'name="MassHunter Data Acquisition" value=""/>\n'
                '    </software>\n'
                '    <software id="ingest" version="1.0">\n'
                '      <cvParam cvRef="MS" accession="MS:1000799" '
                'name="custom unreleased software tool" '
                'value="agilent_d.py"/>\n'
                '    </software>\n  </softwareList>\n'
                % _esc(self.ds.contents.get("AcqSoftwareVersion", "unknown")))

        self._w('  <instrumentConfigurationList count="1">\n'
                '    <instrumentConfiguration id="IC1">\n'
                '      <cvParam cvRef="MS" accession="MS:1000490" '
                'name="Agilent instrument model" value=""/>\n'
                '      <userParam name="model number" value="%s"/>\n'
                '      <userParam name="serial number" value="%s"/>\n'
                '      <componentList count="3">\n'
                '        <source order="1">\n'
                '          <cvParam cvRef="MS" accession="MS:1000389" '
                'name="electron ionization" value=""/>\n'
                '        </source>\n'
                '        <analyzer order="2">\n'
                '          <cvParam cvRef="MS" accession="MS:1000081" '
                'name="quadrupole" value=""/>\n'
                '        </analyzer>\n'
                '        <detector order="3">\n'
                '          <cvParam cvRef="MS" accession="MS:1000253" '
                'name="electron multiplier" value=""/>\n'
                '        </detector>\n'
                '      </componentList>\n'
                '      <softwareRef ref="acquisition"/>\n'
                '    </instrumentConfiguration>\n'
                '  </instrumentConfigurationList>\n'
                % (_esc(model), _esc(serial)))

        self._w('  <dataProcessingList count="1">\n'
                '    <dataProcessing id="DP1">\n'
                '      <processingMethod order="1" softwareRef="ingest">\n'
                '        <cvParam cvRef="MS" accession="MS:1000544" '
                'name="Conversion to mzML" value=""/>\n'
                '      </processingMethod>\n'
                '    </dataProcessing>\n  </dataProcessingList>\n')

        started = self.ds.acquired_time or ""
        self._w('  <run id="%s" defaultInstrumentConfigurationRef="IC1" '
                'sampleRef="SAMPLE1" startTimeStamp="%s">\n'
                % (run_id, _esc(started)))

    def _spectrum_list(self):
        ds = self.ds
        n = len(ds.scans)
        self._w('    <spectrumList count="%d" '
                'defaultDataProcessingRef="DP1">\n' % n)
        for i, r in enumerate(ds.scans):
            native = "scanId=%d" % r.scan_id
            self._w('      ')
            # mzML index offsets must point at the '<spectrum' byte itself.
            self.spec_offsets.append((native, self._pos()))
            mz, ab = ds.spectrum(i)
            self._w('<spectrum index="%d" id="%s" defaultArrayLength="%d">\n'
                    % (i, native, len(mz)))
            self._w('        <cvParam cvRef="MS" accession="MS:1000579" '
                    'name="MS1 spectrum" value=""/>\n')
            self._w('        <cvParam cvRef="MS" accession="MS:1000511" '
                    'name="ms level" value="%d"/>\n' % max(1, r.ms_level))
            self._w('        <cvParam cvRef="MS" accession="MS:1000127" '
                    'name="centroid spectrum" value=""/>\n')
            self._w('        <cvParam cvRef="MS" accession="MS:1000%s" '
                    'name="%s scan" value=""/>\n'
                    % ("130" if r.ion_polarity == 0 else "129",
                       "positive" if r.ion_polarity == 0 else "negative"))
            self._w('        <cvParam cvRef="MS" accession="MS:1000285" '
                    'name="total ion current" value="%s"/>\n' % _num(r.tic))
            if mz:
                self._w('        <cvParam cvRef="MS" accession="MS:1000504" '
                        'name="base peak m/z" value="%s" unitCvRef="MS" '
                        'unitAccession="MS:1000040" unitName="m/z"/>\n'
                        % _num(r.base_peak_mz))
                self._w('        <cvParam cvRef="MS" accession="MS:1000505" '
                        'name="base peak intensity" value="%s" '
                        'unitCvRef="MS" unitAccession="MS:1000131" '
                        'unitName="number of detector counts"/>\n'
                        % _num(r.base_peak_value))
                self._w('        <cvParam cvRef="MS" accession="MS:1000528" '
                        'name="lowest observed m/z" value="%s" '
                        'unitCvRef="MS" unitAccession="MS:1000040" '
                        'unitName="m/z"/>\n' % _num(r.min_x))
                self._w('        <cvParam cvRef="MS" accession="MS:1000527" '
                        'name="highest observed m/z" value="%s" '
                        'unitCvRef="MS" unitAccession="MS:1000040" '
                        'unitName="m/z"/>\n' % _num(r.max_x))
            self._w('        <scanList count="1">\n'
                    '          <cvParam cvRef="MS" accession="MS:1000795" '
                    'name="no combination" value=""/>\n'
                    '          <scan>\n'
                    '            <cvParam cvRef="MS" accession="MS:1000016" '
                    'name="scan start time" value="%s" unitCvRef="UO" '
                    'unitAccession="UO:0000031" unitName="minute"/>\n'
                    % _num(r.scan_time))
            self._w('            <scanWindowList count="1">\n'
                    '              <scanWindow>\n'
                    '                <cvParam cvRef="MS" '
                    'accession="MS:1000501" name="scan window lower limit" '
                    'value="%s" unitCvRef="MS" unitAccession="MS:1000040" '
                    'unitName="m/z"/>\n'
                    '                <cvParam cvRef="MS" '
                    'accession="MS:1000500" name="scan window upper limit" '
                    'value="%s" unitCvRef="MS" unitAccession="MS:1000040" '
                    'unitName="m/z"/>\n'
                    '              </scanWindow>\n'
                    '            </scanWindowList>\n'
                    '          </scan>\n        </scanList>\n'
                    % (_num(r.mass_range_min), _num(r.mass_range_max)))
            self._w('        <binaryDataArrayList count="2">\n')
            self._binary_array(mz, "mz", "          ")
            self._binary_array(ab, "intensity", "          ")
            self._w('        </binaryDataArrayList>\n')
            self._w('      </spectrum>\n')
        self._w('    </spectrumList>\n')

    def _chromatogram_list(self):
        rt, sig = self.ds.tic()
        self._w('    <chromatogramList count="1" '
                'defaultDataProcessingRef="DP1">\n')
        self._w('      ')
        self.chrom_offsets.append(("TIC", self._pos()))
        self._w('<chromatogram index="0" id="TIC" '
                'defaultArrayLength="%d">\n' % len(rt))
        self._w('        <cvParam cvRef="MS" accession="MS:1000235" '
                'name="total ion current chromatogram" value=""/>\n')
        self._w('        <binaryDataArrayList count="2">\n')
        self._binary_array(rt, "time", "          ")
        self._binary_array(sig, "intensity", "          ")
        self._w('        </binaryDataArrayList>\n')
        self._w('      </chromatogram>\n    </chromatogramList>\n')

    def _index(self):
        self._w('<indexList count="2">\n')
        self._w('  <index name="spectrum">\n')
        for native, off in self.spec_offsets:
            self._w('    <offset idRef="%s">%d</offset>\n' % (_esc(native), off))
        self._w('  </index>\n  <index name="chromatogram">\n')
        for native, off in self.chrom_offsets:
            self._w('    <offset idRef="%s">%d</offset>\n' % (_esc(native), off))
        self._w('  </index>\n</indexList>\n')

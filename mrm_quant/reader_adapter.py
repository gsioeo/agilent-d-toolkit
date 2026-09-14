"""Read-only adapter over the existing ``ingest/agilent_d.py`` reader.

All binary access stays in the existing reader: this module loads it by explicit
file path, so nothing depends on the working directory, checks that the dataset
really has the layout that reader assumes, and converts its scan records into
the canonical dictionaries the rest of the package consumes.

Canonical record keys
---------------------
``scan_id``, ``rt_min``, ``scan_type``, ``ms_level``, ``precursor_mz``,
``collision_energy_ev``, ``polarity``, ``time_segment_id``, ``scan_method_id``,
``cycle_number``, ``product_mz``, ``intensity``, ``declared_product_mz``, plus
``tic`` and ``point_count`` carried through for diagnostics.

Intensities are the values stored by the instrument. They are not converted to
counts or counts per second, and the dwell time is not applied; the unit is
reported everywhere as ``stored_intensity``.
"""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET

#: Spectrum format identifiers seen in the validated datasets. A format ID is
#: never treated as sufficient on its own; the byte counts, offsets and point
#: counts are checked as well.
DEFAULT_ALLOWED_FORMAT_IDS = (1, 3)

#: Scan types that carry a selected precursor, i.e. that can hold a transition.
MRM_SCAN_TYPES = (256,)

#: XSD primitive -> stored width in bytes, used only to confirm that the record
#: stride the existing reader assumes matches the schema shipped with the data.
_XSD_WIDTHS = {'int': 4, 'long': 8, 'double': 8, 'float': 4,
               'short': 2, 'byte': 1, 'unsignedInt': 4, 'unsignedShort': 2}

_READER_CACHE = {}


class UnsupportedLayout(ValueError):
    """Raised when a dataset does not match the validated binary layout."""


def default_ingest_dir():
    """The ``ingest/`` directory of this repository, or ``MRM_QUANT_INGEST_DIR``."""
    env = os.environ.get('MRM_QUANT_INGEST_DIR')
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parents[1] / 'ingest'


def load_reader(ingest_dir=None):
    """Import ``agilent_d.py`` read-only, by explicit path, and cache it."""
    directory = Path(ingest_dir) if ingest_dir else default_ingest_dir()
    source = (directory / 'agilent_d.py').resolve()
    cached = _READER_CACHE.get(str(source))
    if cached is not None:
        return cached
    if not source.is_file():
        raise ValueError('reader_not_found: %s' % source)
    spec = importlib.util.spec_from_file_location('mrm_quant_legacy_agilent_d', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _READER_CACHE[str(source)] = module
    return module


def _local(tag):
    return tag.rsplit('}', 1)[-1]


def _schema_record_size(acqdata):
    """Record stride declared by ``MSScan.xsd``, with one SpectrumParams group."""
    path = Path(acqdata) / 'MSScan.xsd'
    if not path.is_file():
        raise UnsupportedLayout('unsupported_layout: MSScan.xsd is missing, '
                                'the record layout cannot be confirmed')
    raw = path.read_bytes()
    if raw.startswith(b'\xef\xbb\xbf'):
        raw = raw[3:]
    root = ET.fromstring(raw)
    types = {}
    for node in root:
        if _local(node.tag) == 'complexType' and node.get('name'):
            types[node.get('name')] = node

    def width(node, seen):
        total = 0
        for element in node.iter():
            if _local(element.tag) != 'element' or not element.get('type'):
                continue
            kind = element.get('type').rsplit(':', 1)[-1]
            if kind in _XSD_WIDTHS:
                total += _XSD_WIDTHS[kind]
            elif kind in types:
                if kind in seen:
                    raise UnsupportedLayout('unsupported_layout: recursive type %s' % kind)
                total += width(types[kind], seen | {kind})
            else:
                raise UnsupportedLayout('unsupported_layout: unknown field type %s' % kind)
        return total

    if 'ScanRecordType' not in types:
        raise UnsupportedLayout('unsupported_layout: MSScan.xsd declares no ScanRecordType')
    return width(types['ScanRecordType'], {'ScanRecordType'})


def _spectrum_params_groups(acqdata):
    """Number of SpectrumParamValues elements declared per scan record."""
    path = Path(acqdata) / 'MSScan.xsd'
    raw = path.read_bytes()
    if raw.startswith(b'\xef\xbb\xbf'):
        raw = raw[3:]
    root = ET.fromstring(raw)
    count = 0
    for element in root.iter():
        if _local(element.tag) == 'element' and element.get('name') == 'SpectrumParamValues':
            count += 1
    return count


def validate_layout(dataset, reader, *, allowed_format_ids=None):
    """Confirm the dataset matches the validated layout; raise otherwise.

    Returns a dictionary of the facts that were checked, for provenance.
    """
    allowed = tuple(allowed_format_ids if allowed_format_ids is not None
                    else DEFAULT_ALLOWED_FORMAT_IDS)
    acqdata = Path(dataset.acqdata)
    peak_path = acqdata / 'MSPeak.bin'
    if not peak_path.is_file():
        raise UnsupportedLayout('unsupported_layout: MSPeak.bin is missing; '
                                'profile-only storage is not supported')
    stride = getattr(reader, '_REC').size
    schema_size = _schema_record_size(acqdata)
    if schema_size != stride:
        raise UnsupportedLayout(
            'unsupported_layout: MSScan.xsd implies a %d byte record, the reader '
            'assumes %d' % (schema_size, stride))
    groups = _spectrum_params_groups(acqdata)
    if groups != 1:
        raise UnsupportedLayout(
            'unsupported_layout: %d SpectrumParamValues groups per record' % groups)

    scans = dataset.scans
    if not scans:
        raise UnsupportedLayout('unsupported_layout: no scan records')
    scan_size = (acqdata / 'MSScan.bin').stat().st_size
    if (scan_size - len(scans) * stride) < 0:
        raise UnsupportedLayout('unsupported_layout: MSScan.bin shorter than its record array')
    peak_size = peak_path.stat().st_size
    end = 0
    for record in scans:
        if record.point_count < 0 or record.byte_count != record.point_count * 8:
            raise UnsupportedLayout(
                'unsupported_layout: scan %d declares %d bytes for %d points'
                % (record.scan_id, record.byte_count, record.point_count))
        if record.spectrum_offset < 0 or record.spectrum_offset + record.byte_count > peak_size:
            raise UnsupportedLayout(
                'unsupported_layout: scan %d peak block runs outside MSPeak.bin'
                % record.scan_id)
        if record.spectrum_format_id not in allowed:
            raise UnsupportedLayout(
                'unsupported_layout: scan %d uses spectrum format %d'
                % (record.scan_id, record.spectrum_format_id))
        end = max(end, record.spectrum_offset + record.byte_count)
    if end != peak_size:
        raise UnsupportedLayout(
            'unsupported_layout: last peak block ends at %d, MSPeak.bin is %d bytes'
            % (end, peak_size))

    declared = sum(int(segment.get('NumOfScans') or 0)
                   for segment in dataset.time_segments)
    if declared and declared != len(scans):
        raise UnsupportedLayout(
            'unsupported_layout: MSTS.xml declares %d scans, MSScan.bin holds %d'
            % (declared, len(scans)))
    return {'record_stride_bytes': stride, 'schema_record_bytes': schema_size,
            'spectrum_params_groups': groups, 'scan_count': len(scans),
            'peak_file_bytes': peak_size, 'scan_file_bytes': scan_size,
            'spectrum_format_ids': sorted({r.spectrum_format_id for r in scans}),
            'allowed_format_ids': list(allowed)}


def channel_id(channel):
    """Stable, filesystem-safe identifier for a method channel."""
    if channel['channel_type'] == 'mrm':
        return 'seg%s_m%s_%gto%g_ce%g' % (
            channel['time_segment_id'], channel['scan_method_id'],
            channel['precursor_mz'], channel['product_mz'],
            channel['collision_energy_ev'])
    return 'seg%s_m%s_ms1_%gto%g' % (
        channel['time_segment_id'], channel['scan_method_id'],
        channel['mz_low'], channel['mz_high'])


def _float(value):
    if value is None or value == '':
        return None
    return float(value)


def _channels_from_method(dataset):
    """Method-declared channels, one row per MRM transition or MS1 segment."""
    polarity_by_method = {}
    for record in dataset.scans:
        key = (record.time_segment_id, record.scan_method_id)
        polarity_by_method.setdefault(key, record.ion_polarity)
    scans_by_method = {}
    for record in dataset.scans:
        key = (record.time_segment_id, record.scan_method_id)
        scans_by_method[key] = scans_by_method.get(key, 0) + 1

    channels = []
    for segment in dataset.ms_method.get('timeSegments', []):
        segment_id = int(segment.get('index') or 0)
        for scan_segment in segment.get('scanSegments', []):
            method_id = int(scan_segment.get('index') or 0)
            key = (segment_id, method_id)
            scan_type = (scan_segment.get('scanType') or '').strip()
            common = dict(time_segment_id=segment_id, scan_method_id=method_id,
                          method_scan_type=scan_type,
                          data_storage=scan_segment.get('dataStorage'),
                          polarity=polarity_by_method.get(key),
                          frame_count=scans_by_method.get(key, 0))
            elements = scan_segment.get('scanElements', [])
            if scan_type.upper() == 'MRM':
                for element in elements:
                    channel = dict(common)
                    channel.update(
                        channel_type='mrm',
                        element_index=int(element.get('index') or 0),
                        compound_name=element.get('compoundName'),
                        precursor_mz=_float(element.get('ms1LowMz')),
                        product_mz=_float(element.get('ms2LowMz')),
                        precursor_resolution=element.get('ms1Res'),
                        product_resolution=element.get('ms2Res'),
                        collision_energy_ev=_float(element.get('collisionEnergy')),
                        dwell_ms=_float(element.get('dwell')),
                        gain=_float(element.get('gain')),
                        is_istd=(element.get('isISTD') == 'true'))
                    channel['channel_id'] = channel_id(channel)
                    channels.append(channel)
            else:
                for element in elements or [{}]:
                    channel = dict(common)
                    channel.update(
                        channel_type='ms1',
                        element_index=int(element.get('index') or 0),
                        mz_low=_float(element.get('ms1LowMz')),
                        mz_high=_float(element.get('ms1HighMz')),
                        step_size=_float(element.get('ms1Stepsize')),
                        gain=_float(element.get('gain')),
                        scan_time_ms=_float(scan_segment.get('scanTime')))
                    channel['channel_id'] = channel_id(channel)
                    channels.append(channel)
    return channels


def _declared_products(channels):
    """(time_segment_id, scan_method_id) -> ascending declared product m/z."""
    declared = {}
    for channel in channels:
        if channel['channel_type'] != 'mrm':
            continue
        key = (channel['time_segment_id'], channel['scan_method_id'])
        declared.setdefault(key, []).append(channel['product_mz'])
    return {key: sorted(value) for key, value in declared.items()}


class RunReader:
    """One dataset, opened read-only, exposing canonical records and method info."""

    def __init__(self, path, *, ingest_dir=None, allowed_format_ids=None, validate=True):
        self.path = Path(path)
        self.reader = load_reader(ingest_dir)
        self.dataset = self.reader.AgilentDotD(str(self.path))
        self.layout = (validate_layout(self.dataset, self.reader,
                                       allowed_format_ids=allowed_format_ids)
                       if validate else {})
        self.channels = _channels_from_method(self.dataset)
        self._declared = _declared_products(self.channels)

    # -- lifecycle ---------------------------------------------------------
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def close(self):
        self.dataset.close()

    # -- metadata ----------------------------------------------------------
    @property
    def run_name(self):
        return self.dataset.name

    def metadata(self):
        dataset = self.dataset
        scans = dataset.scans
        return {'run_name': dataset.name, 'dataset_path': str(self.path),
                'sample_name': dataset.sample_name,
                'acquired_time': dataset.acquired_time,
                'instrument': dataset.instrument,
                'method_name': dataset.method_name,
                'scan_count': len(scans),
                'mrm_frame_count': sum(1 for s in scans if s.scan_type in MRM_SCAN_TYPES),
                'ms1_frame_count': sum(1 for s in scans if s.ms_level == 1),
                'rt_min_first': scans[0].scan_time if scans else None,
                'rt_min_last': scans[-1].scan_time if scans else None,
                'method_fingerprint': self.method_fingerprint(),
                'intensity_unit': 'stored_intensity'}

    def method_signature(self):
        """The method facts a batch must share, in a stable, comparable form."""
        method = self.dataset.ms_method
        channels = []
        for channel in sorted(self.channels, key=lambda c: (c['time_segment_id'],
                                                            c['scan_method_id'],
                                                            c['element_index'])):
            channels.append({key: channel[key] for key in sorted(channel)
                             if key not in ('frame_count',)})
        return {'ms_instrument': method.get('msInstrument'),
                'ion_source': method.get('ionSource'),
                'stop_mode': method.get('stopMode'),
                'stop_time': method.get('stopTime'),
                'solvent_delay': method.get('solventDelay'),
                'collision_gas_on': method.get('collisionGasOn'),
                'use_gain': method.get('useGain'),
                'time_filter': method.get('isTimeFilterEnabled'),
                'time_filter_peak_width': method.get('timeFilterPeakWidth'),
                'source_parameters': method.get('sourceParameters'),
                'time_segments': [{'index': s.get('index'),
                                   'start_time': s.get('startTime')}
                                  for s in method.get('timeSegments', [])],
                'channels': channels}

    def method_fingerprint(self):
        """Short digest of :meth:`method_signature`.

        Derived from the parsed acquisition method, not from the method file
        name and not from the printed ``acqmeth.txt``, so two runs match only
        when the acquired channels and source settings agree.
        """
        blob = json.dumps(self.method_signature(), sort_keys=True,
                          ensure_ascii=False, default=str).encode('utf-8')
        return 'M-' + hashlib.sha256(blob).hexdigest()[:16]

    # -- records -----------------------------------------------------------
    def header(self, record):
        return {'scan_id': record.scan_id, 'rt_min': record.scan_time,
                'scan_type': record.scan_type, 'ms_level': record.ms_level,
                'precursor_mz': record.mz_of_interest,
                'collision_energy_ev': record.collision_energy,
                'polarity': record.ion_polarity,
                'time_segment_id': record.time_segment_id,
                'scan_method_id': record.scan_method_id,
                'cycle_number': record.cycle_number,
                'tic': record.tic, 'point_count': record.point_count}

    def iter_records(self, header_filter=None):
        """Canonical records, decoding peak lists only for the kept frames."""
        for index, record in enumerate(self.dataset.scans):
            head = self.header(record)
            if header_filter is not None and not header_filter(head):
                continue
            mz, intensity = self.dataset.spectrum(index)
            head['product_mz'] = [float(value) for value in mz]
            head['intensity'] = [float(value) for value in intensity]
            head['declared_product_mz'] = list(self._declared.get(
                (record.time_segment_id, record.scan_method_id), []))
            yield head

    def records(self, header_filter=None):
        return list(self.iter_records(header_filter))


def read_records(path, *, header_filter=None, **options):
    with RunReader(path, **options) as run:
        return run.records(header_filter)


def method_channels(path, **options):
    with RunReader(path, **options) as run:
        return run.channels


def method_fingerprint(path, **options):
    with RunReader(path, **options) as run:
        return run.method_fingerprint()


def run_metadata(path, **options):
    with RunReader(path, **options) as run:
        return run.metadata()

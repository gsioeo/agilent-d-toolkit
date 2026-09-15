# Agilent MassHunter `.d` binary layout (GC/MS, MSScan.bin + MSPeak.bin)

Reverse-engineered and validated against every scan of a nine-run GC/MS
sequence. Written for anyone who needs to read Agilent GC/MS data without
MassHunter or the vendor .NET data-access DLLs.

Layout observed on data written by MassHunter Acquisition 13.x from a
triple-quadrupole GC/MS in single-quad (MS1 full scan) mode, centroid storage.
It should be checked, not assumed, on other instrument and software
combinations; the consistency tests at the end of this document are how to do
that.

Reference implementation: `ingest/agilent_d.py`.

## Dataset tree

```
<run>.d/
  acqmeth.txt                     human-readable method printout, UTF-16LE
  sequence.log                    cumulative instrument sequence log
  runstart.txt  GC.INI  ...       GC state at run start
  AcqData/
    Contents.xml                  acquisition time, software version, technique
    Devices.xml                   MS model / serial / firmware
    sample_info.xml               vial, injection volume, tune file, operator
    sequence.xml                  the sequence row that produced this run
    MSTS.xml                      time segments: start/end RT, scan count
    MSScan.xsd                    field order of the MSScan.bin records
    MSScan.bin                    fixed-size scan index records
    MSPeak.bin                    centroid peak lists
    <method>.m/                   full method copy incl. qqqacqmethod.xml
    *.pdf                         tune and method reports
```

`MSScan.xsd` is the key: it declares the exact field sequence of a scan record.
The binary is that sequence packed little-endian with no alignment padding.

## MSScan.bin

```
offset 0x00   first bytes are 01 01 00 00 (format marker)
offset 0x58   int32    byte offset where the record array starts (296 here)
offset <that> record array, 184 bytes per scan, count = (filesize - start) / 184
```

Record layout (all little-endian, no padding). Offsets are within the record:

| off | type   | field                      | example value           |
|-----|--------|----------------------------|-------------------------|
| 0   | int32  | ScanID                     | 1                       |
| 4   | int32  | ScanMethodID               | 1                       |
| 8   | int32  | TimeSegmentID              | 1                       |
| 12  | double | ScanTime (minutes)         | 10.003183333333334      |
| 20  | int32  | MSLevel                    | 1                       |
| 24  | int32  | ScanType                   | 1 (Scan)                |
| 28  | double | TIC                        | 165111.79136657715      |
| 36  | double | BasePeakMZ                 | 595.0999755859375       |
| 44  | double | BasePeakValue              | 5156.4619140625         |
| 52  | int32  | CycleNumber                | 1                       |
| 56  | int32  | Status                     | 0                       |
| 60  | int32  | IonMode                    | 2 (EI)                  |
| 64  | int32  | IonPolarity                | 0 (positive)            |
| 68  | double | Fragmentor                 | 0.0                     |
| 76  | double | CollisionEnergy            | 0.0                     |
| 84  | double | MzOfInterest               | 0.0                     |
| 92  | double | SamplingPeriod             | 0.0                     |
| 100 | double | MeasuredMassRangeMin       | 50.0                    |
| 108 | double | MeasuredMassRangeMax       | 600.0                   |
| 116 | double | Threshold                  | 100.0                   |
| 124 | int32  | IsFragmentorDynamic        | 0                       |
| 128 | int32  | IsCollisionEnergyDynamic   | 0                       |
| 132 | int32  | SpectrumFormatID           | 1                       |
| 136 | int64  | SpectrumOffset             | 68                      |
| 144 | int32  | ByteCount                  | 816                     |
| 148 | int32  | PointCount                 | 102                     |
| 152 | double | MinY (lowest abundance)    | 104.13257598876953      |
| 160 | double | MaxY (highest abundance)   | 5156.4619140625         |
| 168 | double | MinX (lowest m/z)          | 506.70001220703125      |
| 176 | double | MaxX (highest m/z)         | 598.5                   |

Total 184 bytes. The struct format string is
`"<3i d 2i 3d 4i 7d 2i i q 2i 4d"`.

Fields 0..128 are `ScanRecordType` from the XSD; fields 132..176 are one
`SpectrumParamsType`. A record carries exactly one spectrum parameter block in
this data. If a future dataset carries more, the record size will not be 184
and the consistency check below will catch it.

`IonMode` and `ScanType` are bit flags; the labels used by the reader are in
`agilent_d.SCAN_TYPES` / `ION_MODES` / `ION_POLARITY`.

## MSPeak.bin

```
offset 0x00   first bytes are 03 01 00 00 (format marker)
offset 0x44   first peak block (matches SpectrumOffset of scan 1)
```

For each scan, at `SpectrumOffset`:

```
PointCount x float32    m/z, strictly ascending
PointCount x float32    abundance
```

so `ByteCount == PointCount * 8`. Blocks are stored in scan order and pack
tightly: the last block ends exactly at EOF.

These are **centroids**, not profile data. This method used
`dataStorage = PeakDetected` with `threshold = 100`, so only centroids at or
above 100 counts were kept. That is why `MinX`/`MaxX` are much narrower than
the 50-600 scan range in early scans, and why a scan may hold ~100 points at
the solvent-delay boundary and ~500 mid-run.

m/z values sit on the 0.1 u acquisition grid (`ms1Stepsize`), stored as
float32. Values such as 82.9 and 96.9 are real: the centroid positions are not
integers.

## Reading it

```python
import sys; sys.path.insert(0, "ingest")
from agilent_d import AgilentDotD

with AgilentDotD("data/run1.d") as ds:
    print(ds.sample_name, ds.instrument, ds.method_name, len(ds.scans))
    rt, tic = ds.tic()                       # per-scan retention time and TIC
    mz, ab = ds.spectrum(0)                  # centroids of the first scan
    rt, traces = ds.eic([93, 107], tol=0.3)  # extracted ion chromatograms
    ds.write_mzml("run1.mzML")               # indexed mzML 1.1.0
```

`AgilentDotD` is standard library only. `eic()` / `flat_peaks()` need numpy;
`xic(mz_low, mz_high)` is the pure-python equivalent.

## Validation

Every check below was run over every scan in the private validation datasets,
not on a sample:

| check | result |
|-------|--------|
| `sum(abundances) == TIC` stored in the record | exact, every scan |
| `mz[argmax(ab)] == BasePeakMZ`, `max(ab) == BasePeakValue` | exact, every scan |
| m/z array strictly ascending | holds, every scan |
| `mz[0] == MinX`, `mz[-1] == MaxX` | exact, every scan |
| `ByteCount == PointCount * 8` | holds, every scan |
| last peak block ends at MSPeak.bin EOF | exact, every file |
| record count vs `NumOfScans` in MSTS.xml | matches, every file |
| first/last `ScanTime` vs MSTS.xml `StartTime`/`EndTime` | matches |

The TIC identity is the strong one: it is an independent quantity stored in a
different file from the peak lists, so reproducing it across the full validation set
confirms the record stride, the peak-block offsets, the float32 interpretation
and the m/z-then-abundance ordering all at once.

## Pitfalls

**Text encoding.** `acqmeth.txt`, `sequence.log` and `runstart.txt` are
UTF-16LE with a BOM. Reading them as UTF-8 or a single-byte codec produces
mojibake; on this instrument the labels are Chinese.
`agilent_d.decode_text_file()` sniffs the BOM.

**Do not assume the record array starts at 296.** Read the int32 at 0x58.
`agilent_d._probe_data_offset()` is a fallback that locates the array by
consistency check if that value is ever wrong.

**Sample name is not the folder name.** `sample_info.xml` and `sequence.log`
can disagree with the directory, and `sequence.log` can list a third set of
names again, because that log is a cumulative record of the whole instrument
sequence rather than of one output folder. Decide explicitly which of the three
is authoritative for your purpose; this toolkit uses the folder name as the
label and records the `sample_info.xml` name alongside it.

**`sample_info.xml` keys.** Some `<Field>` entries have a numeric hash in
`<Name>`; the only human-readable key is `<DisplayName>`, which is localised
(on a non-English instrument it is not in English). The reader falls back to
`DisplayName` when `Name` is all digits.

**m/z window edges.** With a 0.1 u grid stored as float32, a request like
`81 +/- 0.3` puts its edges exactly on the 80.7 and 81.3 grid points, where
float rounding decides inclusion. Worse, NumPy's NEP-50 promotion casts a
float64 bound down to float32 when compared against a float32 array, so a
vectorised comparison silently includes points a pure-python one excludes.
The reader widens both mass arrays to float64 and adds `EIC_EDGE = 1e-4` to
each side so edge grid points are always included, deterministically.

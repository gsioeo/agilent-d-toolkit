# `ingest/` toolkit reference

Four tools plus a driver, reading Agilent MassHunter `.d` datasets directly.
Binary layout is documented in [agilent-d-format.md](agilent-d-format.md).
MRM transition extraction and same-batch calibration quantification live in a
separate package, documented in [mrm-quant.md](mrm-quant.md); it reuses this
reader read-only and does not change anything below.

## Install

Nothing to install if the files are present. To set the toolkit up in another
dataset folder:

```bash
cd /path/to/new/dataset/folder
patch -p1 < gcms-agilent-toolkit.patch
python3 ingest/run_all.py --mzml --verify
```

Requirements: python3 >= 3.8. `agilent_d.py` is standard library only; numpy is
needed for `eic.py` and `plot.py`, matplotlib for the figures.

## Files

| file | role |
|------|------|
| `ingest/agilent_d.py` | reader library: metadata, scans, spectra, EIC, mzML writer |
| `ingest/ingest.py` | stage 1: manifest, per-scan TIC, summed spectra, mzML |
| `ingest/plot.py` | stage 2: TIC and summed-spectrum figures |
| `ingest/eic.py` | stage 3: extracted ion chromatograms, peak table, figures |
| `ingest/run_all.py` | runs all three stages with one set of options |

## One command

```bash
python3 ingest/run_all.py                 # metadata + figures + EIC, no mzML
python3 ingest/run_all.py --mzml          # add mzML export (~45 MB per run)
python3 ingest/run_all.py --only run1 run2   # a subset of runs
python3 ingest/run_all.py --mz 93 107 --tol 0.3
python3 ingest/run_all.py --skip eic      # stop after the figures
python3 ingest/run_all.py --source /other/data --out /tmp/out --verify
```

`--verify` runs the EIC extraction self-checks (see Validation).

## Output tree

```
ingested/
  manifest.json                    40 metadata fields per run
  samples.csv                      the same, one row per run
  method/<method>.m/               verbatim copy of the acquisition method
  method/<method>.acqmeth.txt      method printout decoded UTF-16LE -> UTF-8
  data/<run>/tic.csv               per scan: RT, TIC, base peak, peak count, m/z range
  data/<run>/avg_spectrum.csv      summed spectrum, 0.1 u bins
  data/<run>/avg_spectrum_unit.csv summed spectrum, 1 u bins
  mzml/<run>.mzML                  indexed mzML 1.1.0, all scans + TIC chromatogram
  plots/tic_grid.png               all runs, common intensity scale
  plots/tic_overlay.png            overlay plus normalised stack
  plots/tic/<run>.png              linear and log panels, peak RTs labelled
  plots/avg_spectrum/<run>.png     full range plus 24 u zoom on the base peak
  plots/avg_spectrum_unit/<run>.png  full range linear and log
  plots/avg_spectrum_unit_grid.png all runs at 1 u
  eic/<run>.csv                    RT, TIC, one eic_<mz> column per ion
  eic/peaks.csv                    run, mz, rt, height, area, width, boundaries
  eic/plots/<run>.png              TIC plus every ion, stacked
  eic/plots/mz<mz>.png             one ion across all runs, common scale
```

Rough sizes per run: CSV ~0.6 MB, figures ~0.4 MB, mzML ~45 MB.

## Stage options worth knowing

`ingest.py`
- `--mzml` / `--no-compress` — mzML export, zlib-compressed float32 by default.
- `--only NAME ...` — restrict runs; matches with or without the `.d` suffix.

`plot.py`
- `--logy` — log intensity on the multi-run spectrum grid. Per-run figures
  always carry both a linear and a second panel.
- `--dpi`, `--only`.

`eic.py`
- `--mz 71 57 85` — explicit targets. Without it, `--auto N` (default 6) picks
  the N most intense unit-mass ions across the selected runs, so the choice
  depends on which runs are included.
- `--tol` — half-window in u, default 0.3. See the m/z edge note below.
- `--rt LO HI` — restrict output and figures to a retention-time window.
- `--max-peaks`, `--no-plots`, `--verify`.

## Library use

```python
import sys; sys.path.insert(0, "ingest")
from agilent_d import AgilentDotD, find_datasets

for path in find_datasets("."):
    with AgilentDotD(path) as ds:
        print(ds.sample_name, ds.acquired_time, len(ds.scans))
        rt, tic = ds.tic()
        mz, ab = ds.spectrum(0)
        rt, traces = ds.eic([71, 93], tol=0.3)     # needs numpy
        rows = ds.average_spectrum(decimals=1, rt_range=(14.4, 14.9))
        ds.write_mzml("/tmp/%s.mzML" % ds.name)
```

Useful members: `sample_info`, `contents`, `devices`, `time_segments`,
`ms_method`, `method_text()`, `scans` (list of `ScanRecord`), `spectrum(i)`,
`iter_spectra()`, `tic()`, `xic(lo, hi)`, `eic(targets, tol)`, `eic_window()`,
`flat_peaks()`, `average_spectrum()`, `write_mzml()`.

## Validation

Run `--verify` to reproduce the EIC checks. The following were confirmed on
every scan of every run of a nine-run validation sequence.

Reader, every scan of every run:
- decoded abundances sum to the stored `TIC` exactly
- abundance argmax reproduces `BasePeakMZ` and `BasePeakValue` exactly
- m/z arrays strictly ascending and equal to `MinX`/`MaxX` at the ends
- last peak block ends exactly at `MSPeak.bin` EOF
- scan count matches `MSTS.xml`

mzML, every file:
- validates against the official `mzML1.1.2_idx.xsd`
- re-read with pymzML and compared back to the raw binaries: max |Δm/z| = 0,
  max |ΔI| = 0, max |ΔRT| = 0
- every index offset lands on a `<spectrum` byte; `fileChecksum` matches the
  SHA-1 of the prefix

EIC:
- a single window spanning 50-600 reproduces the per-scan TIC with max
  difference exactly 0
- vectorised and pure-python extraction agree bit-for-bit
- no EIC point exceeds its scan's TIC

Portable patch: applies cleanly into an empty folder, restores all five sources
byte-identically, and the full pipeline then runs with self-checks passing.

## Gotchas

**m/z window edges.** The mass axis is a 0.1 u grid stored as float32, so
`81 +/- 0.3` lands its edges exactly on grid points. NumPy's NEP-50 promotion
also casts a float64 bound down to float32 when compared against a float32
array, which silently includes points a pure-python comparison excludes. The
reader widens the arrays to float64 and adds `EIC_EDGE = 1e-4` per side. Ask
for `--tol 0.3` and you get `[t-0.3001, t+0.3001]`, i.e. grid points at exactly
`+/-0.3` are included.

**Unit windows do not tile the grid.** Extracting all 551 unit masses at
`--tol 0.45` leaves ~1.8e-4 of the signal out: nothing covers the `x.5` grid
points. Verified that the residual equals the abundance on those points
exactly. Use one wide window, not a sum of narrow ones, if you need closure.

**mzML ids must not start with a digit.** A run whose folder name begins with a
digit yields an `id` beginning with a digit, which fails the mzML NCName
pattern. Names are sanitised
(`_ncname()`), and the original sample name is carried in `sampleList`.

**Auto-picked EIC ions depend on the run selection**, because they come from
the summed spectra of the selected runs. Pass `--mz` for a fixed ion set when
comparing across sessions.

**Reported values in mzML cvParams round-trip exactly.** `%g`-style formatting
was losing precision (20 microseconds on scan start time); values are written
with `repr()`.

**`ingested/` is regenerated through a temporary stage.** A completed stage
replaces only artifacts it generated, preserving unexpected user files in the
output directory; an interrupted stage leaves the prior output intact.
Duplicate `.d` basenames are rejected before any processing. `manifest.json`
only lists an `mzml` output if that run was ingested with `--mzml` in the same
invocation.

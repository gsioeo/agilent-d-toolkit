# agilent-d-toolkit

**Read Agilent MassHunter GC/MS `.d` data without MassHunter.** A pure-Python
reader for `MSScan.bin` / `MSPeak.bin`, with indexed mzML export, TIC and
extracted-ion-chromatogram extraction, publication-ready plots, and MRM
transition extraction with same-batch calibration quantification.

`Python 3.8+` · `no third-party deps for the reader` · `mzML 1.1.0`

---

## Why

Agilent `.d` is an undocumented binary format. The usual ways to read it are
MassHunter itself, the vendor .NET data-access DLLs, or a Windows COM bridge —
all of which mean a licence, a specific OS, or both. That is a poor fit for a
Linux analysis box, a container, or a reproducible pipeline.

This toolkit reads the binaries directly. The layout was derived from the
`MSScan.xsd` schema that Agilent ships inside every dataset and then confirmed
numerically against the instrument's own stored values, so it is not guesswork —
see [Correctness](#correctness).

## Install

No build step. Either copy `ingest/` next to your data, or carry the single
patch file across and unpack it in place:

```bash
cp gcms-agilent-toolkit.patch /path/to/folder/containing/your/.d/directories
cd /path/to/folder/containing/your/.d/directories
patch -p1 < gcms-agilent-toolkit.patch
```

The reader is standard library only; numpy and matplotlib are needed for the EIC
stage and the figures.

## Quick start

```bash
python3 ingest/run_all.py --mzml --verify
```

That walks every `.d` directory below the current folder and writes everything
to `ingested/`. On a 9-run, 240 MB sequence the full pipeline including mzML
export takes about 50 seconds.

```bash
python3 ingest/run_all.py                     # skip the mzML export
python3 ingest/run_all.py --only run1 run2    # a subset of runs
python3 ingest/run_all.py --mz 93 107 --tol 0.3
python3 ingest/run_all.py --skip eic          # stop after the figures
python3 ingest/run_all.py --source /data --out /tmp/out
```

`--verify` runs the extraction self-checks on your own data and should report a
maximum difference of exactly 0 against the TIC stored in the raw files.

## R plotting

After the Python ingestion has written `ingested/`, you can make ggplot2-based
figures without re-reading the raw Agilent binaries:

```bash
Rscript ingest/plot_r.R
Rscript ingest/plot_r.R --only STD s22 --dpi 300
```

The R script reads the same `ingested/data/<run>/` and optional
`ingested/eic/` CSV files produced by `python3 ingest/run_all.py`, then writes
PNG figures under `ingested/plots_r/`. It uses base R plus `ggplot2`; if
`ggplot2` is missing, the script installs it from CRAN before plotting.

## MRM quantification

`mrm_quant/` is a separate package for MRM data: it extracts one acquired
transition, integrates it inside an expected retention-time window, fits a
calibration line from the standards of the **same batch** and back-calculates
the samples. It reuses the reader above read-only and is not wired into
`run_all.py`.

```bash
python3 -m mrm_quant inspect  --source ../batch-folder --out quant_results/inspection_001
python3 -m mrm_quant extract   --batch batch.csv --analytes analytes.json --out quant_results/traces_001
python3 -m mrm_quant quantify  --batch batch.csv --analytes analytes.json --out quant_results/quant_001
```

`inspect` needs no concentrations: it lists every acquired channel with its
Q1/Q3, collision energy, dwell and frame count, and ranks the channels by
response so the quantifier can be chosen. That choice then has to be written
into `analytes.json` — `quantify` never picks a channel by response, and never
falls back to an MS1 extracted-ion chromatogram at the same nominal mass. A
transition the method did not acquire is an error (`transition_not_acquired`).

Two UTF-8 configuration files drive it; templates are in
`mrm_quant/templates/`. `batch.csv` gives one row per injection and analyte with
an explicit `role` (calibration / blank / qc / unknown), the level, the
concentration and the dilution factor. `analytes.json` gives the transition, the
retention-time window, the integration settings, the calibration model
(`none`, `1/x` or `1/x2` weighting; free or zero intercept) and the QC limits.
Concentrations may be written level by level, generated from a series
(`top_concentration`, `dilution_step`, `levels`), or both — and then they must
agree.

```
quant_results/quant_001/
  integration.csv           per run: selected peak, bounds, baseline, area, response
  calibration_points.csv    per level: concentration, response, back-calculation, bias
  calibration_models.json   slope, intercept, equation, r2, r2_weighted, range
  results.csv               status, vial and reported concentration, unit, validated
  qc.csv                    qualifier ratios, blank check, independent QC
  plots/                    calibration line and one peak figure per run
  provenance.json           raw SHA-256, config hash, method fingerprint, versions
```

Only a clean result reports a concentration. Everything else keeps its
diagnostic numbers and a status instead: `no_peak`, `ambiguous_peak`,
`non_positive_area`, `negative_backcalc`, `below_calibration_range`,
`above_calibration_range`, `below_validated_loq`, `ion_ratio_fail`,
`rt_mismatch`, `blank_contamination`, `qc_fail`. Nothing is extrapolated,
clipped to zero or reported as ND. `validated` is true only when a blank limit
and an independent QC were configured and passed; an absent check is reported as
not evaluated, never as a pass. A model belongs to one batch and one method
fingerprint, so another batch cannot borrow it.

Intensities are the values the instrument stored (`stored_intensity`); the dwell
time is not applied and no vendor comparison has been made, so areas
(`stored_intensity*min`) are comparable within a batch rather than absolute.
[docs/mrm-quant.md](docs/mrm-quant.md) has the full configuration reference and
the limits.

## What you get

```
ingested/
  manifest.json                     40 metadata fields per run
  samples.csv                       the same, one row per run
  method/<method>.acqmeth.txt       method printout, UTF-16LE decoded to UTF-8
  method/<method>.m/                verbatim copy of the acquisition method
  data/<run>/tic.csv                per scan: RT, TIC, base peak, peak count
  data/<run>/avg_spectrum.csv       summed spectrum, 0.1 u bins
  data/<run>/avg_spectrum_unit.csv  summed spectrum, 1 u bins
  mzml/<run>.mzML                   indexed mzML 1.1.0, all scans + TIC
  plots/...                         TIC and spectrum figures, per run and pooled
  eic/<run>.csv                     RT, TIC, one column per target ion
  eic/peaks.csv                     RT, height, trapezoidal area, width, bounds
  eic/plots/...                     EIC per run and per ion across runs
```

Roughly 2 MB of CSV and figures per run, plus ~45 MB per run if you export mzML.

## Library use

```python
import sys; sys.path.insert(0, "ingest")
from agilent_d import AgilentDotD, find_datasets

for path in find_datasets("."):
    with AgilentDotD(path) as ds:
        print(ds.name, ds.acquired_time, len(ds.scans))

        rt, tic = ds.tic()                          # per-scan RT and TIC
        mz, ab  = ds.spectrum(0)                    # centroids of one scan
        rt, tr  = ds.eic([71, 93], tol=0.3)         # ion chromatograms (numpy)
        rows    = ds.average_spectrum(decimals=1, rt_range=(14.4, 14.9))

        ds.write_mzml("/tmp/%s.mzML" % ds.name)
```

Scan records are named tuples, so `ds.scans[i].scan_time`, `.tic`,
`.base_peak_mz`, `.point_count` and the rest are all directly available.
`xic(mz_low, mz_high)` is the numpy-free equivalent of `eic()`.

## Correctness

`MSScan.bin` stores a TIC per scan. `MSPeak.bin` stores the centroids. They are
independent, so reconstructing the stored TIC from the decoded peak lists tests
the record stride, the block offsets, the float32 interpretation and the
array ordering simultaneously. On the validation sequence this matched
**exactly for every scan of every run** — as did the base-peak m/z and
intensity, the m/z ordering, and the observed mass range.

Also verified:

- every generated mzML validates against the official `mzML1.1.2_idx.xsd`
- re-reading with pymzML and comparing back to the raw binaries gives a maximum
  absolute difference of **zero** in m/z, intensity and retention time
- all index offsets resolve to a `<spectrum` element and the `fileChecksum`
  matches the SHA-1 of the preceding bytes
- one EIC window spanning the full mass range reproduces the per-scan TIC with a
  difference of exactly zero
- the vectorised and pure-python extraction paths agree bit for bit

[docs/toolkit.md](docs/toolkit.md) has the full table, plus the traps worth
knowing: m/z window edges landing on the 0.1 u grid, NumPy NEP-50 scalar
promotion silently narrowing a window, and the mzML id constraints.
[docs/agilent-d-format.md](docs/agilent-d-format.md) documents the byte layout
field by field.

`mrm_quant/` is specified by 28 executable contracts written before the code and
run against the real data: the September batch's `204 -> 93` channel must return
its 2454 native MRM points with the stored first-frame intensity, MS1 frames may
not leak into a transition, a declared product missing from a frame is `None`
rather than zero, a measured zero stays zero, and a mass tolerance wide enough to
cover two acquired channels is rejected instead of merged. The integrator is
checked against hand-computed areas on unevenly spaced points, and the three
weighting options against hand-computed slopes and intercepts; the reported `r2`
was compared with `numpy.polyfit` plus Pearson r² on the real curve and agrees to
the last digit. The whole suite is 115 tests:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

## Scope and limits

The reader and the `ingest/` pipeline were developed and validated against
**MS1 full-scan EI data, centroid (`PeakDetected`) storage, one time segment**,
written by MassHunter Acquisition 13.x from a triple-quadrupole GC/MS.
`mrm_quant/` adds **MRM (`ScanType` 256, MS level 2)** on two batches from the
same instrument, where MRM and MS1 frames alternate within one cycle.

Not exercised, and likely to need work:

- **SIM acquisitions.** The scan record carries `MzOfInterest`,
  `CollisionEnergy` and `ScanType`, and they are parsed, but no SIM dataset was
  available to test against.
- **Scheduled MRM with several time segments or several collision energies per
  transition.** The selector keys exist and are matched, and two acquisition
  groups are refused rather than summed, but only single-segment data was seen.
- **Profile data.** `MSProfile.bin` is not read at all; only centroids.
- **Multiple time segments.** Scans from every segment are read and exported,
  but the manifest cross-checks its scan count against the first segment only.
- **LC/MS `.d` variants.** Untouched. The record layout comes from the same
  schema, so it may work, but the mzML writer hardcodes EI as the ion source.
- **Vendor agreement.** No MassHunter export was compared against, so areas and
  concentrations are internally consistent rather than verified against the
  vendor's own integrator.

Out of scope by design: no library search or compound identification, no
deconvolution, no baseline correction beyond a straight line between the
integration bounds. The `ingest/` peak areas are plain trapezoids between the
local minima either side of the apex (cut off at 5% of peak height), which is
fine for comparing runs and not a substitute for a real integrator;
`mrm_quant/` integrates the same way but inside a configured window, with
explicit bounds, interpolated endpoints and a status for every peak it will not
report.

## Requirements

| | |
|---|---|
| Python | >= 3.8 |
| numpy | optional — needed by `eic.py` and `plot.py` |
| matplotlib | optional — needed for the figures |
| R | optional — needed by `ingest/plot_r.R` |
| ggplot2 | optional — installed by `plot_r.R` if missing |

`agilent_d.py` on its own has no third-party dependencies.

## Repository layout

```
ingest/agilent_d.py   reader library and mzML writer
ingest/ingest.py      stage 1: manifest, TIC, summed spectra, mzML
ingest/plot.py        stage 2: TIC and summed-spectrum figures
ingest/plot_r.R       optional R/ggplot2 plotting from ingested CSVs
ingest/eic.py         stage 3: extracted ion chromatograms, peak table
ingest/run_all.py     all three stages under one set of options
ingest/make_patch.py  regenerates gcms-agilent-toolkit.patch
mrm_quant/            MRM extraction and same-batch quantification (v1)
  reader_adapter.py   read-only adapter over ingest/agilent_d.py, layout checks
  transitions.py      transition selection, MS1 TIC and MRM sum kept apart
  integration.py      windowed integration and retention-time peak selection
  calibration.py      weighted linear fit, response, back-calculation
  qc.py               qualifier ratio, blank, status aggregation, output guard
  config.py           batch.csv and analytes.json, concentration series
  pipeline.py         inspect, extract, quantify
  report.py           tables, provenance, figures
  templates/          example configuration, no sample identities
tests/                115 tests: legacy regression, 28 contracts, real batches
docs/                 format reference, tool reference, mrm-quant reference
```

`gcms-agilent-toolkit.patch` is a self-contained installer holding every
`ingest/` source and its docs as a plain unified diff, so a new dataset folder
needs no clone — just `patch -p1`. Run `python3 ingest/make_patch.py` after
editing anything there so the two cannot drift. `mrm_quant/` ships with the
repository and is not part of that patch.

## Data policy

This repository contains **no instrument data**. `.gitignore` excludes `.d`
directories, raw instrument files, tune and method reports, generated output,
figures, CSV and mzML, plus the quantification inputs and outputs
(`quant_config/`, `batch*.csv`, `analytes*.json`, `quant_results/`) because a
batch table names samples and concentrations. Instrument serial numbers,
operator names, sample identities and acquisition paths exist only in the raw
files, in `ingested/` and in those local configurations, none of which is
tracked. Before publishing, check `git status --ignored`.

> **Careful:** if the repository lives in the same folder as your data, that
> data is *ignored*, not *tracked* — `git clean -fdx` would delete it. Use
> `git clean -fd` without `-x`.

## Licence

[GNU Affero General Public License v3.0](LICENSE) (`AGPL-3.0-only`).

Copyleft, including over a network: if you run a modified version as a service,
the AGPL requires you to offer that version's source to its users. If you intend
this code to be embeddable in closed pipelines, a more permissive licence would
suit better.

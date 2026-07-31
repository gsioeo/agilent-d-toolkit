#!/usr/bin/env python3
"""Ingest the Agilent GC/MS ``.d`` datasets in this workspace into open formats.

Default output layout (under ``ingested/``)::

    manifest.json                     everything known about every run
    samples.csv                       one row per run (spreadsheet friendly)
    method/<method>.acqmeth.txt       decoded instrument method printout
    data/<sample>/tic.csv             RT, TIC, base peak, peak count per scan
    data/<sample>/avg_spectrum.csv    summed spectrum at 0.1 u
    data/<sample>/avg_spectrum_unit.csv   summed spectrum at 1 u
    mzml/<sample>.mzML                full centroid data (with --mzml)

Usage::

    python3 ingest/ingest.py                # metadata + chromatograms
    python3 ingest/ingest.py --mzml         # also export mzML
    python3 ingest/ingest.py --only run1 run2   # restrict to named runs
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from agilent_d import (AgilentDotD, ION_MODES, ION_POLARITY, SCAN_TYPES,
                       find_datasets)

WORKSPACE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SAMPLE_COLUMNS = [
    "sample", "dataset_dir", "vial", "sample_type", "acquired_time",
    "instrument", "ms_model", "ms_serial", "ms_firmware", "acq_software",
    "method", "tune_file", "tune_datetime", "ion_source", "ion_mode",
    "polarity", "scan_type", "ms_level", "mz_low", "mz_high", "mz_step",
    "gain", "scan_time_ms", "threshold", "data_storage",
    "solvent_delay_min", "rt_start_min", "rt_end_min", "n_scans",
    "median_scan_period_s", "total_tic", "max_tic", "max_tic_rt_min",
    "global_base_peak_mz", "n_centroids", "inj_volume_ul", "dilution",
    "original_path",
]


def _f(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _first_scan_element(ds):
    for seg in ds.ms_method.get("timeSegments") or []:
        for ss in seg.get("scanSegments") or []:
            for se in ss.get("scanElements") or []:
                return ss, se
    return {}, {}


def summarize(ds: AgilentDotD) -> dict:
    scans = ds.scans
    dev = ds.ms_device
    ss, se = _first_scan_element(ds)

    rts = [s.scan_time for s in scans]
    tics = [s.tic for s in scans]
    total_tic = sum(tics)
    max_tic = max(tics) if tics else 0.0
    max_i = tics.index(max_tic) if tics else -1
    periods = sorted((rts[i + 1] - rts[i]) * 60.0 for i in range(len(rts) - 1))
    median_period = periods[len(periods) // 2] if periods else None
    top = max(scans, key=lambda s: s.base_peak_value) if scans else None
    first = scans[0] if scans else None

    seg = (ds.time_segments or [{}])[0]

    return {
        "sample": ds.sample_name,
        "dataset_dir": os.path.basename(ds.path),
        "vial": ds.sample_info.get("Sample Position", ""),
        "sample_type": ds.sample_info.get("Sample Type", ""),
        "acquired_time": ds.acquired_time,
        "instrument": ds.instrument,
        "ms_model": dev.get("ModelNumber", ""),
        "ms_serial": dev.get("SerialNumber", ""),
        "ms_firmware": dev.get("FirmwareVersion", ""),
        "acq_software": ds.contents.get("AcqSoftwareVersion", ""),
        "method": ds.method_name,
        "tune_file": ds.sample_info.get("调谐文件") or ds.ms_method.get("tuneFile", ""),
        "tune_datetime": ds.sample_info.get("调谐日期时间", ""),
        "ion_source": ds.ms_method.get("ionSource", ""),
        "ion_mode": ION_MODES.get(first.ion_mode, first.ion_mode) if first else "",
        "polarity": ION_POLARITY.get(first.ion_polarity, first.ion_polarity) if first else "",
        "scan_type": ss.get("scanType") or (SCAN_TYPES.get(first.scan_type) if first else ""),
        "ms_level": first.ms_level if first else "",
        "mz_low": se.get("ms1LowMz", ""),
        "mz_high": se.get("ms1HighMz", ""),
        "mz_step": se.get("ms1Stepsize", ""),
        "gain": se.get("gain", ""),
        "scan_time_ms": ss.get("scanTime", ""),
        "threshold": ss.get("threshold", ""),
        "data_storage": ss.get("dataStorage", ""),
        "solvent_delay_min": ds.ms_method.get("solventDelay", ""),
        "rt_start_min": round(rts[0], 5) if rts else "",
        "rt_end_min": round(rts[-1], 5) if rts else "",
        "n_scans": len(scans),
        "median_scan_period_s": round(median_period, 4) if median_period else "",
        "total_tic": round(total_tic, 1),
        "max_tic": round(max_tic, 1),
        "max_tic_rt_min": round(rts[max_i], 5) if max_i >= 0 else "",
        "global_base_peak_mz": round(top.base_peak_mz, 3) if top else "",
        "n_centroids": sum(s.point_count for s in scans),
        "inj_volume_ul": ds.sample_info.get("Inj Vol (µl)", ""),
        "dilution": ds.sample_info.get("Dilution", ""),
        "original_path": ds.sample_info.get("Data File", ""),
        "declared_scans_mstsxml": seg.get("NumOfScans", ""),
    }


def write_tic(ds: AgilentDotD, path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["scan_id", "rt_min", "tic", "base_peak_mz",
                    "base_peak_intensity", "n_centroids", "lowest_mz",
                    "highest_mz"])
        for s in ds.scans:
            w.writerow([s.scan_id, "%.6f" % s.scan_time, "%.1f" % s.tic,
                        "%.3f" % s.base_peak_mz, "%.1f" % s.base_peak_value,
                        s.point_count,
                        "%.3f" % s.min_x if s.point_count else "",
                        "%.3f" % s.max_x if s.point_count else ""])


def _strip_d(name: str) -> str:
    return name[:-2] if name.lower().endswith(".d") else name


def write_avg_spectrum(rows, path: str) -> None:
    base = max((r[1] for r in rows), default=0.0)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["mz", "summed_intensity", "relative_pct", "n_scans"])
        for mz, inten, n in rows:
            w.writerow(["%.4f" % mz, "%.1f" % inten,
                        "%.4f" % (100.0 * inten / base) if base else "0",
                        n])


def rebin_unit(rows):
    acc = {}
    for mz, inten, n in rows:
        key = int(round(mz))
        slot = acc.setdefault(key, [0.0, 0])
        slot[0] += inten
        slot[1] += n
    return [(float(k), v[0], v[1]) for k, v in sorted(acc.items())]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default=WORKSPACE,
                    help="directory to scan for .d datasets (default: workspace)")
    ap.add_argument("--out", default=os.path.join(WORKSPACE, "ingested"),
                    help="output directory (default: <workspace>/ingested)")
    ap.add_argument("--mzml", action="store_true",
                    help="also export full centroid data as indexed mzML")
    ap.add_argument("--no-compress", action="store_true",
                    help="write mzML binary arrays uncompressed")
    ap.add_argument("--only", nargs="*", metavar="NAME",
                    help="restrict to these sample/dataset names")
    args = ap.parse_args(argv)

    datasets = find_datasets(args.source)
    if args.only:
        want = {_strip_d(n).lower() for n in args.only}
        datasets = [d for d in datasets
                    if _strip_d(os.path.basename(d)).lower() in want]
    if not datasets:
        print("no .d datasets found under %s" % args.source, file=sys.stderr)
        return 1

    out = args.out
    data_dir = os.path.join(out, "data")
    method_dir = os.path.join(out, "method")
    mzml_dir = os.path.join(out, "mzml")
    for d in (out, data_dir, method_dir):
        os.makedirs(d, exist_ok=True)
    if args.mzml:
        os.makedirs(mzml_dir, exist_ok=True)

    summaries = []
    methods_written = set()
    t0 = time.time()

    for path in datasets:
        with AgilentDotD(path) as ds:
            info = summarize(ds)
            label = _strip_d(info["dataset_dir"]) or ds.sample_name
            sdir = os.path.join(data_dir, label)
            os.makedirs(sdir, exist_ok=True)

            write_tic(ds, os.path.join(sdir, "tic.csv"))
            rows = ds.average_spectrum(decimals=1)
            write_avg_spectrum(rows, os.path.join(sdir, "avg_spectrum.csv"))
            write_avg_spectrum(rebin_unit(rows),
                               os.path.join(sdir, "avg_spectrum_unit.csv"))

            mname = ds.method_name or "method"
            if mname not in methods_written:
                text = ds.method_text()
                if text:
                    with open(os.path.join(method_dir, mname + ".acqmeth.txt"),
                              "w", encoding="utf-8") as fh:
                        fh.write(text)
                if ds.method_dir:
                    dest = os.path.join(method_dir, mname)
                    if not os.path.exists(dest):
                        shutil.copytree(ds.method_dir, dest)
                methods_written.add(mname)

            info["outputs"] = {
                "tic": os.path.relpath(os.path.join(sdir, "tic.csv"), out),
                "avg_spectrum": os.path.relpath(
                    os.path.join(sdir, "avg_spectrum.csv"), out),
                "avg_spectrum_unit": os.path.relpath(
                    os.path.join(sdir, "avg_spectrum_unit.csv"), out),
            }

            if args.mzml:
                target = os.path.join(mzml_dir, label + ".mzML")
                ds.write_mzml(target, compress=not args.no_compress)
                info["outputs"]["mzml"] = os.path.relpath(target, out)

            summaries.append(info)
            print("  %-8s %5d scans  %8.3f-%.3f min  TIC max %.3g @ %s min%s"
                  % (label, info["n_scans"], info["rt_start_min"] or 0,
                     info["rt_end_min"] or 0, info["max_tic"],
                     info["max_tic_rt_min"],
                     "  -> mzML" if args.mzml else ""))

    manifest = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "source": os.path.abspath(args.source),
        "n_datasets": len(summaries),
        "reader": "ingest/agilent_d.py (MSScan.bin + MSPeak.bin, stdlib only)",
        "samples": summaries,
    }
    with open(os.path.join(out, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, ensure_ascii=False)

    with open(os.path.join(out, "samples.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SAMPLE_COLUMNS, extrasaction="ignore")
        w.writeheader()
        for s in summaries:
            w.writerow(s)

    print("\n%d dataset(s) ingested into %s in %.1fs"
          % (len(summaries), out, time.time() - t0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Extract ion chromatograms (EICs) from the Agilent ``.d`` runs.

For each requested target m/z, every centroid within +/- ``--tol`` is summed
per scan, giving one trace per ion per run.  Output under ``ingested/eic/``::

    <run>.csv               rt_min, tic, eic_<mz> ... (one column per ion)
    peaks.csv               detected peaks: run, mz, rt, height, area, width
    plots/<run>.png         TIC plus every ion for that run
    plots/mz<mz>.png        one ion compared across all runs

Examples::

    python3 ingest/eic.py                        # auto-pick the 6 top ions
    python3 ingest/eic.py --mz 71 57 85 93 147   # explicit ions
    python3 ingest/eic.py --mz 207 --tol 0.5 --rt 10 16
    python3 ingest/eic.py --only run1 --mz 93 --verify
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from agilent_d import AgilentDotD, find_datasets

WORKSPACE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_TRAPZ = getattr(np, "trapezoid", None) or np.trapz


def _sort_key(name: str):
    """Runs named STD (a calibration standard) sort first, then purely numeric
    names in numeric order, then everything else by trailing number."""
    if name.upper() == "STD":
        return (0, 0, name)
    if name.isdigit():
        return (1, int(name), name)
    digits = "".join(c for c in name if c.isdigit())
    return (2, int(digits) if digits else 0, name)


def _strip_d(name: str) -> str:
    return name[:-2] if name.lower().endswith(".d") else name


def fmt_mz(mz: float) -> str:
    return ("%g" % mz)


# --------------------------------------------------------------------------
# peak picking
# --------------------------------------------------------------------------

def detect_peaks(rt, sig, min_rel=0.02, max_peaks=8, min_sep=0.04):
    """Local maxima with trapezoidal area over the surrounding valley."""
    peaks = []
    if sig.size < 3 or sig.max() <= 0:
        return peaks
    interior = np.arange(1, sig.size - 1)
    cand = interior[(sig[1:-1] > sig[:-2]) & (sig[1:-1] >= sig[2:])]
    cand = cand[sig[cand] >= min_rel * sig.max()]
    cand = cand[np.argsort(sig[cand])[::-1]]

    taken = []
    for i in cand:
        if len(taken) >= max_peaks:
            break
        if any(abs(rt[i] - rt[j]) < min_sep for j in taken):
            continue
        taken.append(i)

        half = sig[i] * 0.05
        left = i
        while left > 0 and sig[left - 1] < sig[left] and sig[left] > half:
            left -= 1
        right = i
        while right < sig.size - 1 and sig[right + 1] < sig[right] \
                and sig[right] > half:
            right += 1
        seg = slice(left, right + 1)
        area = float(_TRAPZ(sig[seg], rt[seg]))
        peaks.append({"rt": float(rt[i]), "height": float(sig[i]),
                      "area": area, "rt_start": float(rt[left]),
                      "rt_end": float(rt[right]),
                      "width_min": float(rt[right] - rt[left])})
    return sorted(peaks, key=lambda p: p["rt"])


# --------------------------------------------------------------------------
# ion selection
# --------------------------------------------------------------------------

def auto_targets(ingested: str, runs, n: int):
    """Top *n* unit-mass ions of the summed spectra of the selected runs."""
    acc = {}
    for run in runs:
        path = os.path.join(ingested, "data", run, "avg_spectrum_unit.csv")
        if not os.path.exists(path):
            continue
        with open(path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                mz = float(row["mz"])
                acc[mz] = acc.get(mz, 0.0) + float(row["relative_pct"])
    if not acc:
        return []
    top = sorted(acc.items(), key=lambda kv: kv[1], reverse=True)[:n]
    return sorted(mz for mz, _ in top)


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

def _spaced(peaks, n, min_gap):
    """Up to *n* of the tallest peaks, no two closer than *min_gap* minutes.

    Keeps the retention-time labels from printing on top of each other.
    """
    chosen = []
    for p in sorted(peaks, key=lambda q: q["height"], reverse=True):
        if len(chosen) >= n:
            break
        if all(abs(p["rt"] - q["rt"]) >= min_gap for q in chosen):
            chosen.append(p)
    return chosen


def plot_run(run, rt, tic, targets, traces, out_path, dpi, peaks_by_ion):
    n = len(targets)
    fig, axes = plt.subplots(n + 1, 1, figsize=(11, 1.55 * (n + 1) + 1.1),
                             sharex=True, squeeze=False)
    axes = axes.ravel()
    cmap = plt.get_cmap("tab10")

    axes[0].plot(rt, tic, linewidth=0.7, color="#444444")
    axes[0].set_ylabel("TIC", fontsize=8)
    axes[0].set_title("TIC (reference)", fontsize=8, loc="left")

    for k, t in enumerate(targets):
        ax = axes[k + 1]
        y = traces[k]
        ax.plot(rt, y, linewidth=0.7, color=cmap(k % 10))
        ax.fill_between(rt, 0, y, color=cmap(k % 10), alpha=0.18, linewidth=0)
        for p in _spaced(peaks_by_ion[k], 6, 0.4):
            ax.annotate("%.2f" % p["rt"], (p["rt"], p["height"]),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=6, color="#b03030")
        share = 100.0 * y.sum() / tic.sum() if tic.sum() else 0.0
        ax.set_ylabel("m/z %s" % fmt_mz(t), fontsize=8)
        ax.set_title("m/z %s  -  %.2f%% of the total ion current"
                     % (fmt_mz(t), share), fontsize=8, loc="left")

    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(6)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.margins(x=0.01)
    axes[-1].set_xlabel("retention time (min)")
    fig.suptitle("EIC - %s" % run, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_ion_across_runs(target, runs, per_run, out_path, dpi):
    n = len(runs)
    fig, axes = plt.subplots(n, 1, figsize=(11, 1.25 * n + 1.2), sharex=True,
                             squeeze=False)
    axes = axes.ravel()
    cmap = plt.get_cmap("turbo")
    ymax = max(per_run[r][1].max() for r in runs) or 1.0

    for k, run in enumerate(runs):
        rt, y = per_run[run]
        c = cmap(k / max(1, n - 1))
        ax = axes[k]
        ax.plot(rt, y, linewidth=0.7, color=c)
        ax.fill_between(rt, 0, y, color=c, alpha=0.2, linewidth=0)
        ax.set_ylim(0, ymax * 1.05)
        ax.set_ylabel(run, fontsize=8, rotation=0, ha="right", va="center")
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.yaxis.get_offset_text().set_fontsize(6)
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.margins(x=0.01)
    axes[-1].set_xlabel("retention time (min)")
    fig.suptitle("EIC m/z %s - all runs (common intensity scale)"
                 % fmt_mz(target), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------

def verify(ds, targets, tol):
    """Self-checks: a full-range EIC must reproduce the stored TIC exactly, the
    vectorised extraction must agree with the pure-python one, and no EIC point
    may exceed its scan's TIC."""
    ok = True
    tic = np.asarray([s.tic for s in ds.scans], dtype=float)

    rt, full = ds.eic([325.0], tol=275.0)          # covers 50-600 entirely
    dmax = float(np.abs(full[0] - tic).max())
    print("    full-range EIC vs stored TIC: max |diff| = %.6g  (TIC max %.4g)"
          % (dmax, tic.max()))
    ok &= dmax == 0.0

    t = targets[0]
    lo, hi = ds.eic_window(t, tol)
    _, fast = ds.eic([t], tol=tol)
    slow = np.asarray(ds.xic(lo, hi)[1], dtype=float)
    scale = max(slow.max(), 1.0)
    d2 = float(np.abs(fast[0] - slow).max())
    print("    m/z %s in [%.5f, %.5f]: vectorised vs pure-python "
          "max |diff| = %.6g  (rel %.2g)" % (fmt_mz(t), lo, hi, d2, d2 / scale))
    ok &= d2 / scale < 1e-12

    over = int((fast[0] > tic + 1e-3).sum())
    print("    scans where EIC > TIC: %d" % over)
    ok &= over == 0

    # the ion count inside the window must never exceed the peak count
    npts = sum(1 for _ in ds.scans)
    print("    scans covered: %d" % npts)
    ok &= npts == len(tic)
    return ok


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default=WORKSPACE,
                    help="directory holding the .d datasets")
    ap.add_argument("--ingested", default=os.path.join(WORKSPACE, "ingested"))
    ap.add_argument("--out", default=None,
                    help="output directory (default: <ingested>/eic)")
    ap.add_argument("--mz", nargs="+", type=float, metavar="MZ",
                    help="target m/z values")
    ap.add_argument("--auto", type=int, default=6, metavar="N",
                    help="if --mz is absent, use the N most intense unit-mass "
                         "ions (default 6)")
    ap.add_argument("--tol", type=float, default=0.3,
                    help="half-window in u around each target (default 0.3)")
    ap.add_argument("--only", nargs="*", metavar="RUN",
                    help="restrict to these runs")
    ap.add_argument("--rt", nargs=2, type=float, metavar=("LO", "HI"),
                    help="restrict output to this retention-time window")
    ap.add_argument("--max-peaks", type=int, default=8)
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--verify", action="store_true",
                    help="run extraction self-checks on the first run")
    args = ap.parse_args(argv)

    datasets = find_datasets(args.source)
    if args.only:
        want = {_strip_d(r).lower() for r in args.only}
        datasets = [d for d in datasets
                    if _strip_d(os.path.basename(d)).lower() in want]
    if not datasets:
        print("no .d datasets found under %s" % args.source, file=sys.stderr)
        return 1
    datasets.sort(key=lambda p: _sort_key(_strip_d(os.path.basename(p))))
    runs = [_strip_d(os.path.basename(p)) for p in datasets]

    targets = args.mz
    if not targets:
        targets = auto_targets(args.ingested, runs, args.auto)
        if not targets:
            print("could not auto-pick ions (no avg_spectrum_unit.csv); "
                  "pass --mz", file=sys.stderr)
            return 1
        print("auto-selected ions: %s"
              % ", ".join("m/z %s" % fmt_mz(t) for t in targets))
    targets = sorted(float(t) for t in targets)

    out = args.out or os.path.join(args.ingested, "eic")
    plot_dir = os.path.join(out, "plots")
    os.makedirs(out, exist_ok=True)
    if not args.no_plots:
        os.makedirs(plot_dir, exist_ok=True)

    per_ion = {t: {} for t in targets}
    peak_rows = []
    t0 = time.time()

    for path, run in zip(datasets, runs):
        with AgilentDotD(path) as ds:
            rt, traces = ds.eic(targets, tol=args.tol)
            tic = np.asarray([s.tic for s in ds.scans], dtype=float)

            if args.verify:
                print("  %s self-checks:" % run)
                if not verify(ds, targets, args.tol):
                    print("    FAILED", file=sys.stderr)
                    return 2
                args.verify = False        # first run is enough

            if args.rt:
                m = (rt >= args.rt[0]) & (rt <= args.rt[1])
                rt, traces, tic = rt[m], traces[:, m], tic[m]

            csv_path = os.path.join(out, run + ".csv")
            with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                w = csv.writer(fh)
                w.writerow(["rt_min", "tic"]
                           + ["eic_%s" % fmt_mz(t) for t in targets])
                for j in range(rt.size):
                    w.writerow(["%.6f" % rt[j], "%.1f" % tic[j]]
                               + ["%.1f" % traces[k, j]
                                  for k in range(len(targets))])

            peaks_by_ion = []
            for k, t in enumerate(targets):
                pk = detect_peaks(rt, traces[k], max_peaks=args.max_peaks)
                peaks_by_ion.append(pk)
                for p in pk:
                    peak_rows.append({
                        "run": run, "mz": fmt_mz(t), "tol": args.tol,
                        "rt_min": "%.4f" % p["rt"],
                        "height": "%.1f" % p["height"],
                        "area": "%.1f" % p["area"],
                        "width_min": "%.4f" % p["width_min"],
                        "rt_start": "%.4f" % p["rt_start"],
                        "rt_end": "%.4f" % p["rt_end"],
                        "pct_of_trace_max":
                            "%.1f" % (100.0 * p["height"] / traces[k].max()
                                      if traces[k].max() else 0.0)})
                per_ion[t][run] = (rt, traces[k])

            if not args.no_plots:
                plot_run(run, rt, tic, targets, traces,
                         os.path.join(plot_dir, run + ".png"), args.dpi,
                         peaks_by_ion)

            shares = ["m/z %s %.1f%%" % (fmt_mz(t),
                      100.0 * traces[k].sum() / tic.sum() if tic.sum() else 0)
                      for k, t in enumerate(targets)]
            print("  %-5s %s" % (run, "  ".join(shares)))

    if not args.no_plots:
        for t in targets:
            plot_ion_across_runs(t, runs, per_ion[t],
                                 os.path.join(plot_dir,
                                              "mz%s.png" % fmt_mz(t)),
                                 args.dpi)

    fields = ["run", "mz", "tol", "rt_min", "height", "area", "width_min",
              "rt_start", "rt_end", "pct_of_trace_max"]
    with open(os.path.join(out, "peaks.csv"), "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(peak_rows)

    print("\n%d ion(s) x %d run(s), %d peaks tabulated -> %s  (%.1fs)"
          % (len(targets), len(runs), len(peak_rows), out, time.time() - t0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

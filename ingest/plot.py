#!/usr/bin/env python3
"""Plot the ingested chromatograms and averaged spectra to PNG.

Reads ``ingested/data/<run>/{tic,avg_spectrum,avg_spectrum_unit}.csv`` and
writes to ``ingested/plots/``::

    tic_grid.png                    all runs, shared axes
    tic_overlay.png                 all runs overlaid, normalised
    tic/<run>.png                   one run, strongest peaks labelled
    avg_spectrum/<run>.png          summed spectrum at 0.1 u
    avg_spectrum_unit/<run>.png     summed spectrum at 1 u
    avg_spectrum_unit_grid.png      all runs at 1 u

Usage::

    python3 ingest/plot.py                 # everything
    python3 ingest/plot.py --only run1 run2   # selected runs
    python3 ingest/plot.py --logy          # log intensity axis on the spectra
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

try:
    from output import begin_output, commit_output, discard_output
except ImportError:
    from ingest.output import begin_output, commit_output, discard_output

WORKSPACE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_IN = os.path.join(WORKSPACE, "ingested")

LINE_KW = dict(linewidth=0.7, color="#1f4e79")


def _sort_key(name: str):
    """Runs named STD (a calibration standard) sort first, then purely numeric
    names in numeric order, then everything else by trailing number."""
    if name.upper() == "STD":
        return (0, 0, name)
    if name.isdigit():
        return (1, int(name), name)
    digits = "".join(c for c in name if c.isdigit())
    return (2, int(digits) if digits else 0, name)


def read_columns(path: str, columns):
    """Read the named float columns of a CSV into numpy arrays."""
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        acc = {c: [] for c in columns}
        for row in reader:
            for c in columns:
                acc[c].append(float(row[c]) if row[c] not in ("", None) else np.nan)
    return {c: np.asarray(v, dtype=float) for c, v in acc.items()}


def find_peaks(x, y, min_rel=0.02, min_sep=0.05, limit=12):
    """Local maxima above *min_rel* of the max, separated by *min_sep* in x."""
    if len(y) < 3:
        return []
    interior = np.arange(1, len(y) - 1)
    cand = interior[(y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:])]
    cand = cand[y[cand] >= min_rel * y.max()]
    cand = cand[np.argsort(y[cand])[::-1]]
    picked = []
    for i in cand:
        if all(abs(x[i] - x[j]) >= min_sep for j in picked):
            picked.append(i)
        if len(picked) >= limit:
            break
    return sorted(picked, key=lambda i: x[i])


def _tidy(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.01)


# --------------------------------------------------------------------------
# TIC
# --------------------------------------------------------------------------

def plot_tic(run, data, out_path, dpi):
    rt, tic = data["rt_min"], data["tic"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.6), sharex=True)

    ax1.plot(rt, tic, **LINE_KW)
    for i in find_peaks(rt, tic, min_rel=0.05, min_sep=0.12, limit=12):
        ax1.annotate("%.2f" % rt[i], (rt[i], tic[i]),
                     textcoords="offset points", xytext=(0, 4), ha="center",
                     fontsize=6.5, color="#b03030")
    ax1.set_ylabel("total ion current")
    ax1.set_ylim(bottom=0)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax1.set_title("linear scale - peaks >5% of base peak labelled",
                  fontsize=9, loc="left")
    _tidy(ax1)

    pos = tic[tic > 0]
    floor = max(pos.min(), np.percentile(pos, 2) / 5.0) if pos.size else 1.0
    ax2.plot(rt, np.clip(tic, floor, None), **LINE_KW)
    ax2.set_yscale("log")
    ax2.set_ylim(floor, tic.max() * 3)
    for i in find_peaks(rt, tic, min_rel=0.004, min_sep=0.3, limit=18):
        ax2.annotate("%.2f" % rt[i], (rt[i], tic[i]),
                     textcoords="offset points", xytext=(0, 4), ha="center",
                     fontsize=6, color="#b03030")
    ax2.set_xlabel("retention time (min)")
    ax2.set_ylabel("total ion current (log)")
    ax2.set_title("log scale - reveals the minor peaks compressed above",
                  fontsize=9, loc="left")
    _tidy(ax2)

    fig.suptitle("TIC - %s  (%d scans, %.2f-%.2f min)"
                 % (run, len(rt), rt[0], rt[-1]), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_tic_grid(runs, tics, out_path, dpi):
    n = len(runs)
    ncols = 3 if n > 4 else max(1, n)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 2.5 * nrows),
                             sharex=True, squeeze=False)
    ymax = max(t["tic"].max() for t in tics.values())
    for ax, run in zip(axes.ravel(), runs):
        d = tics[run]
        ax.plot(d["rt_min"], d["tic"], **LINE_KW)
        ax.set_title(run, fontsize=9)
        ax.set_ylim(0, ymax * 1.05)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.tick_params(labelsize=7)
        ax.yaxis.get_offset_text().set_fontsize(6)
        _tidy(ax)
    for ax in axes.ravel()[len(runs):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("retention time (min)", fontsize=8)
    for row in axes:
        row[0].set_ylabel("TIC", fontsize=8)
    fig.suptitle("TIC - all runs (common intensity scale)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_tic_overlay(runs, tics, out_path, dpi):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    cmap = plt.get_cmap("turbo")
    colors = [cmap(i / max(1, len(runs) - 1)) for i in range(len(runs))]

    for run, c in zip(runs, colors):
        d = tics[run]
        ax1.plot(d["rt_min"], d["tic"], linewidth=0.6, color=c, label=run,
                 alpha=0.85)
    ax1.set_ylabel("total ion current")
    ax1.set_ylim(bottom=0)
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax1.legend(ncol=len(runs), fontsize=7, frameon=False, loc="upper right")
    ax1.set_title("TIC overlay - absolute intensity", fontsize=10)
    _tidy(ax1)

    for k, (run, c) in enumerate(zip(runs, colors)):
        d = tics[run]
        y = d["tic"] / d["tic"].max()
        ax2.plot(d["rt_min"], y + k, linewidth=0.6, color=c)
        ax2.text(d["rt_min"][0] - 0.15, k + 0.35, run, fontsize=7, ha="right",
                 va="center", color=c)
    ax2.set_yticks([])
    ax2.set_xlabel("retention time (min)")
    ax2.set_ylabel("normalised TIC (offset per run)")
    ax2.set_title("TIC stack - each run scaled to its own base peak",
                  fontsize=10)
    ax2.set_xlim(left=tics[runs[0]]["rt_min"][0] - 1.2)
    _tidy(ax2)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------
# spectra
# --------------------------------------------------------------------------

def plot_spectrum(run, data, out_path, dpi, kind, logy=False, n_labels=15):
    mz, rel = data["mz"], data["relative_pct"]
    step = "1 u" if kind == "unit" else "0.1 u"
    min_gap = 7.0 if kind == "unit" else 2.5

    def label(ax, n, gap, max_rel=None, lo=None, hi=None):
        order = np.argsort(rel)[::-1]
        chosen = []
        for i in order:
            if len(chosen) >= n:
                break
            if max_rel is not None and rel[i] >= max_rel:
                continue
            if lo is not None and not (lo <= mz[i] <= hi):
                continue
            if all(abs(mz[i] - mz[j]) >= gap for j in chosen):
                chosen.append(i)
        for i in chosen:
            ax.annotate("%.4g" % mz[i], (mz[i], rel[i]),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=6.5, color="#b03030")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6.6))

    ax1.vlines(mz, 0, rel, linewidth=0.6, color="#1f4e79")
    label(ax1, n_labels, min_gap)
    ax1.set_ylabel("relative intensity (% of base peak)")
    ax1.set_ylim(0, 105)
    ax1.set_xlim(mz.min() - 5, mz.max() + 5)
    ax1.set_xlabel("m/z")
    ax1.set_title("full range, linear - strongest ions labelled", fontsize=9,
                  loc="left")
    _tidy(ax1)

    if kind == "unit":
        # 551 sticks: a log axis is readable and shows the trace-level ions.
        keep = rel > 0.02
        floor = 0.02
        ax2.vlines(mz[keep], floor, rel[keep], linewidth=0.6, color="#1f4e79")
        ax2.set_yscale("log")
        ax2.set_ylim(floor, 400)
        ax2.set_xlim(mz.min() - 5, mz.max() + 5)
        label(ax2, 14, 10.0, max_rel=20.0)
        ax2.set_ylabel("relative intensity (%, log)")
        ax2.set_title("full range, log - strongest ions below 20% labelled "
                      "(bins under 0.02% omitted)", fontsize=9, loc="left")
    else:
        # 5455 sticks would smear into a solid block; zoom instead so the
        # 0.1 u granularity that this file carries is actually visible.
        centre = mz[int(np.argmax(rel))]
        lo, hi = centre - 12.0, centre + 12.0
        keep = (mz >= lo) & (mz <= hi)
        ax2.vlines(mz[keep], 0, rel[keep], linewidth=0.9, color="#1f4e79")
        ax2.set_ylim(0, 105)
        ax2.set_xlim(lo, hi)
        label(ax2, 14, 0.45, lo=lo, hi=hi)
        ax2.set_ylabel("relative intensity (%)")
        ax2.set_title("zoom on the base peak, %.1f-%.1f u - shows the 0.1 u "
                      "bin structure" % (lo, hi), fontsize=9, loc="left")
    ax2.set_xlabel("m/z")
    _tidy(ax2)

    fig.suptitle("Summed spectrum, all scans, %s bins - %s   [%d bins]"
                 % (step, run, len(mz)), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def plot_spectrum_grid(runs, specs, out_path, dpi, logy=False):
    n = len(runs)
    ncols = 3 if n > 4 else max(1, n)
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.6 * ncols, 2.5 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    for ax, run in zip(axes.ravel(), runs):
        d = specs[run]
        ax.vlines(d["mz"], 0, d["relative_pct"], linewidth=0.5,
                  color="#1f4e79")
        top, gap = [], 12.0
        for i in np.argsort(d["relative_pct"])[::-1]:
            if len(top) >= 4:
                break
            if all(abs(d["mz"][i] - d["mz"][j]) >= gap for j in top):
                top.append(i)
        for i in top:
            ax.annotate("%d" % round(d["mz"][i]),
                        (d["mz"][i], d["relative_pct"][i]),
                        textcoords="offset points", xytext=(0, 2),
                        ha="center", fontsize=6, color="#b03030")
        ax.set_title(run, fontsize=9)
        ax.tick_params(labelsize=7)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(1e-3, 200)
        else:
            ax.set_ylim(0, 110)
        _tidy(ax)
    for ax in axes.ravel()[len(runs):]:
        ax.axis("off")
    for ax in axes[-1]:
        ax.set_xlabel("m/z", fontsize=8)
    for row in axes:
        row[0].set_ylabel("rel. int. (%)", fontsize=8)
    fig.suptitle("Summed spectrum, 1 u bins - all runs", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ingested", default=DEFAULT_IN,
                    help="ingested/ directory (default: <workspace>/ingested)")
    ap.add_argument("--out", default=None,
                    help="plot directory (default: <ingested>/plots)")
    ap.add_argument("--only", nargs="*", metavar="RUN",
                    help="restrict to these run names")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--logy", action="store_true",
                    help="log intensity axis on the multi-run spectrum grid "
                         "(per-run figures always show linear + log)")
    args = ap.parse_args(argv)

    data_dir = os.path.join(args.ingested, "data")
    if not os.path.isdir(data_dir):
        print("no ingested data at %s - run ingest/ingest.py first" % data_dir,
              file=sys.stderr)
        return 1

    runs = sorted((d for d in os.listdir(data_dir)
                   if os.path.isdir(os.path.join(data_dir, d))), key=_sort_key)
    if args.only:
        want = {r.lower() for r in args.only}
        runs = [r for r in runs if r.lower() in want]
    if not runs:
        print("no matching runs", file=sys.stderr)
        return 1

    requested_out = args.out or os.path.join(args.ingested, "plots")
    out = begin_output(requested_out)
    try:
        for sub in ("tic", "avg_spectrum", "avg_spectrum_unit"):
            os.makedirs(os.path.join(out, sub), exist_ok=True)

        tics, unit_specs = {}, {}
        written = 0
        for run in runs:
            rd = os.path.join(data_dir, run)

            tic = read_columns(os.path.join(rd, "tic.csv"), ["rt_min", "tic"])
            tics[run] = tic
            plot_tic(run, tic, os.path.join(out, "tic", run + ".png"), args.dpi)

            fine = read_columns(os.path.join(rd, "avg_spectrum.csv"),
                                ["mz", "relative_pct"])
            plot_spectrum(run, fine,
                          os.path.join(out, "avg_spectrum", run + ".png"),
                          args.dpi, kind="fine", logy=args.logy)

            unit = read_columns(os.path.join(rd, "avg_spectrum_unit.csv"),
                                ["mz", "relative_pct"])
            unit_specs[run] = unit
            plot_spectrum(run, unit,
                          os.path.join(out, "avg_spectrum_unit", run + ".png"),
                          args.dpi, kind="unit", logy=args.logy)

            written += 3
            top = np.argsort(unit["relative_pct"])[::-1][:5]
            print("  %-5s TIC max %.3g @ %.2f min | top m/z %s"
                  % (run, tic["tic"].max(),
                     tic["rt_min"][int(np.argmax(tic["tic"]))],
                     ", ".join("%d" % round(unit["mz"][i]) for i in top)))

        plot_tic_grid(runs, tics, os.path.join(out, "tic_grid.png"), args.dpi)
        plot_tic_overlay(runs, tics, os.path.join(out, "tic_overlay.png"), args.dpi)
        plot_spectrum_grid(runs, unit_specs,
                           os.path.join(out, "avg_spectrum_unit_grid.png"),
                           args.dpi, logy=args.logy)
        written += 3

        commit_output(out, requested_out, group='plot')
    except Exception:
        discard_output(out)
        raise

    print("\n%d PNG files written to %s" % (written, requested_out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

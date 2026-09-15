#!/usr/bin/env python3
"""One command to take a folder of Agilent MassHunter ``.d`` runs to results.

Runs the three stages in order:

  1. ingest  - metadata manifest, per-scan TIC, summed spectra, optional mzML
  2. plot    - TIC and summed-spectrum figures
  3. eic     - extracted ion chromatograms, peak table, comparison figures

Everything lands under ``<workspace>/ingested/``.  Stages are independent
scripts (``ingest.py``, ``plot.py``, ``eic.py``); this driver just wires them
together with one set of options.

Usage::

    python3 ingest/run_all.py                     # full pipeline, no mzML
    python3 ingest/run_all.py --mzml              # also export mzML (~45 MB/run)
    python3 ingest/run_all.py --only run1 run2     # a subset of runs
    python3 ingest/run_all.py --mz 93 107 --tol 0.3
    python3 ingest/run_all.py --skip eic          # stop after the plots
    python3 ingest/run_all.py --source /path/to/other/data --out /tmp/out
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import eic as eic_stage
import ingest as ingest_stage
import plot as plot_stage
from agilent_d import find_datasets
try:
    from output import unique_dataset_basenames
except ImportError:
    from ingest.output import unique_dataset_basenames

WORKSPACE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAGES = ("ingest", "plot", "eic")


def _rule(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", default=WORKSPACE,
                    help="directory to scan for .d datasets "
                         "(default: the workspace root)")
    ap.add_argument("--out", default=None,
                    help="output directory (default: <workspace>/ingested)")
    ap.add_argument("--only", nargs="*", metavar="RUN",
                    help="restrict every stage to these run names")
    ap.add_argument("--mzml", action="store_true",
                    help="export full centroid data as indexed mzML")
    ap.add_argument("--mz", nargs="+", type=float, metavar="MZ",
                    help="EIC target m/z (default: 6 auto-picked ions)")
    ap.add_argument("--tol", type=float, default=0.3,
                    help="EIC half-window in u (default 0.3)")
    ap.add_argument("--rt", nargs=2, type=float, metavar=("LO", "HI"),
                    help="restrict the EIC output to this RT window")
    ap.add_argument("--dpi", type=int, default=150)
    ap.add_argument("--skip", nargs="*", default=[], choices=STAGES,
                    metavar="STAGE",
                    help="stages to skip: %s" % ", ".join(STAGES))
    ap.add_argument("--no-plots", action="store_true",
                    help="skip both figure stages (same as --skip plot eic "
                         "for figures, EIC csv still written)")
    ap.add_argument("--verify", action="store_true",
                    help="run the EIC extraction self-checks")
    args = ap.parse_args(argv)

    out = args.out or os.path.join(WORKSPACE, "ingested")
    datasets = find_datasets(args.source)
    unique_dataset_basenames(datasets)
    if not datasets:
        print("no .d datasets found under %s" % args.source, file=sys.stderr)
        return 1

    print("source     : %s" % os.path.abspath(args.source))
    print("output     : %s" % os.path.abspath(out))
    print("datasets   : %d found%s"
          % (len(datasets),
             "" if not args.only else " (restricted to %s)" % " ".join(args.only)))
    print("mzML export: %s" % ("yes" if args.mzml else "no"))

    t0 = time.time()
    common_only = (["--only"] + list(args.only)) if args.only else []

    if "ingest" not in args.skip:
        _rule("stage 1/3  ingest -- metadata, TIC, summed spectra%s"
              % (", mzML" if args.mzml else ""))
        argv1 = ["--source", args.source, "--out", out] + common_only
        if args.mzml:
            argv1.append("--mzml")
        rc = ingest_stage.main(argv1)
        if rc:
            return rc

    if "plot" not in args.skip and not args.no_plots:
        _rule("stage 2/3  plot -- TIC and summed-spectrum figures")
        rc = plot_stage.main(["--ingested", out, "--dpi", str(args.dpi)]
                             + common_only)
        if rc:
            return rc

    if "eic" not in args.skip:
        _rule("stage 3/3  eic -- extracted ion chromatograms")
        argv3 = ["--source", args.source, "--ingested", out,
                 "--tol", str(args.tol), "--dpi", str(args.dpi)] + common_only
        if args.mz:
            argv3 += ["--mz"] + [repr(m) for m in args.mz]
        if args.rt:
            argv3 += ["--rt", str(args.rt[0]), str(args.rt[1])]
        if args.no_plots:
            argv3.append("--no-plots")
        if args.verify:
            argv3.append("--verify")
        rc = eic_stage.main(argv3)
        if rc:
            return rc

    _rule("done in %.1fs" % (time.time() - t0))
    for rel in ("manifest.json", "samples.csv", "data", "method", "mzml",
                "plots", "eic"):
        p = os.path.join(out, rel)
        if os.path.exists(p):
            print("  %s" % p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

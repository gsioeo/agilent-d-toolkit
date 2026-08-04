#!/usr/bin/env python3
"""Regenerate ``gcms-agilent-toolkit.patch`` from the tracked sources.

The patch is a plain ``--- /dev/null`` unified diff, so it needs only
``patch -p1`` to apply -- no git repository at the destination.

Run after changing any file listed in :data:`FILES`::

    python3 ingest/make_patch.py
"""

from __future__ import annotations

import hashlib
import os
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "gcms-agilent-toolkit.patch")

FILES = [
    "ingest/agilent_d.py",
    "ingest/ingest.py",
    "ingest/plot.py",
    "ingest/plot_r.R",
    "ingest/eic.py",
    "ingest/run_all.py",
    "ingest/make_patch.py",
    "docs/agilent-d-format.md",
    "docs/toolkit.md",
]

HEADER = """\
Agilent MassHunter GC/MS .d toolkit -- portable patch
=====================================================

Recreates a self-contained toolkit that reads Agilent MassHunter `.d` datasets
(MSScan.bin + MSPeak.bin) with no vendor software, and produces metadata
manifests, chromatograms, summed spectra, extracted ion chromatograms and
indexed mzML 1.1.0.

Apply into any folder that contains `.d` directories:

    cd /path/to/new/dataset/folder
    patch -p1 < gcms-agilent-toolkit.patch
    python3 ingest/run_all.py --mzml --verify

`--verify` runs the extraction self-checks; they should report a max difference
of exactly 0 against the TIC stored in the raw files.

Contents
    ingest/agilent_d.py       reader library: metadata, scans, spectra, EIC, mzML
    ingest/ingest.py          stage 1: manifest, TIC, summed spectra, mzML
    ingest/plot.py            stage 2: TIC and summed-spectrum figures
    ingest/plot_r.R           optional R/ggplot2 plots from ingested CSVs
    ingest/eic.py             stage 3: extracted ion chromatograms, peak table
    ingest/run_all.py         all three stages under one set of options
    ingest/make_patch.py      regenerates this patch
    docs/agilent-d-format.md  binary layout reference and pitfalls
    docs/toolkit.md           tool reference, output tree, validation, gotchas

Requirements
    python3      >= 3.8       reader is standard library only
    numpy        optional     needed by eic.py and plot.py
    matplotlib   optional     needed for the figures
    R            optional     needed by ingest/plot_r.R
    ggplot2      optional     installed by plot_r.R if missing

Nothing outside the target folder is touched, the source `.d` directories are
only ever opened for reading, and output goes to `ingested/`.

This patch carries code and documentation only -- no instrument data, no sample
information, no local paths.

Generated %(when)s -- %(nfiles)d files, %(nlines)d lines.
%(sums)s
--- 8< --- patch body below, do not edit --- 8< ---

"""


def build() -> str:
    bodies, sums, total = [], [], 0
    for rel in FILES:
        with open(os.path.join(ROOT, rel), "rb") as fh:
            raw = fh.read()
        lines = raw.decode("utf-8").split("\n")
        if lines and lines[-1] == "":
            lines.pop()                 # trailing newline is not a line
        total += len(lines)
        sums.append("    %-26s %5d lines  sha256 %s"
                    % (rel, len(lines), hashlib.sha256(raw).hexdigest()[:32]))
        bodies.append("\n".join(
            ["diff --git a/%s b/%s" % (rel, rel),
             "new file mode 100644",
             "--- /dev/null",
             "+++ b/%s" % rel,
             "@@ -0,0 +1,%d @@" % len(lines)]
            + ["+" + line for line in lines]))

    head = HEADER % {"when": time.strftime("%Y-%m-%d %H:%M:%S%z"),
                     "nfiles": len(FILES), "nlines": total,
                     "sums": "\n".join(sums)}
    return head + "\n".join(bodies) + "\n"


def main() -> int:
    text = build()
    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
    print("wrote %s" % OUT)
    print("  %d bytes, %d files" % (len(text.encode("utf-8")), len(FILES)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

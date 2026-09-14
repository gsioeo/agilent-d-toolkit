"""Path resolution for the repository test suite.

The package under test lives in this repository; the raw ``.d`` datasets live
outside it, under a data root that is configurable through ``GCMS_DATA_ROOT``
and defaults to the repository's parent directory. Nothing here asserts
anything: it only resolves paths and loads the existing reader by explicit
file path, so the tests stay independent of the working directory.
"""
import importlib.util
import os
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
DATA_ROOT = Path(os.environ.get('GCMS_DATA_ROOT') or REPO.parent).resolve()
LEGACY_DIR = REPO / 'ingest'

sys.dont_write_bytecode = True
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def load_legacy(module_name='agilent_d'):
    """Import an existing ``ingest/`` module read-only, by explicit path."""
    spec = importlib.util.spec_from_file_location(
        'legacy_' + module_name, LEGACY_DIR / (module_name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def dataset_path(relative):
    """Resolve a dataset path against the data root, then the repository."""
    for base in (DATA_ROOT, REPO):
        candidate = base / relative
        if candidate.exists():
            return candidate
    return DATA_ROOT / relative

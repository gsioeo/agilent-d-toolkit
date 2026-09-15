"""Shared paths and optional, local-only real-data fixtures for tests."""
import importlib.util
import json
import os
from pathlib import Path
import sys

REPO = Path(__file__).resolve().parents[1]
LEGACY_DIR = REPO / 'ingest'
FIXTURE_MANIFEST = Path(os.environ.get(
    'GCMS_TEST_MANIFEST', REPO / 'tests' / 'data-paths.local.json')).resolve()

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


def fixture_manifest():
    """Load private fixture paths without embedding them in tracked tests."""
    if not FIXTURE_MANIFEST.is_file():
        return {}
    with open(FIXTURE_MANIFEST, encoding='utf-8') as handle:
        document = json.load(handle)
    return document.get('datasets', {})


def dataset_path(key):
    """Return a named private dataset, or a non-existent path when absent."""
    entry = fixture_manifest().get(key, {})
    value = entry.get('path') if isinstance(entry, dict) else entry
    return Path(value).expanduser().resolve() if value else REPO / '.missing-fixture' / key


def fixture_value(key, name, default=None):
    entry = fixture_manifest().get(key, {})
    return entry.get(name, default) if isinstance(entry, dict) else default

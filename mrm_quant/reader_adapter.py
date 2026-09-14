"""Read-only adapter over the existing ingest/agilent_d.py reader.

Loads the legacy reader by explicit file path so nothing depends on the working
directory, and converts its scan records into the canonical dictionaries the
rest of this package consumes.
"""


def load_reader(ingest_dir=None):
    raise NotImplementedError('reader_adapter.load_reader is implemented in task 2')


def read_records(path, ingest_dir=None):
    raise NotImplementedError('reader_adapter.read_records is implemented in task 2')


def method_channels(path, ingest_dir=None):
    raise NotImplementedError('reader_adapter.method_channels is implemented in task 2')

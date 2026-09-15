"""Adapter tests: canonical records, method channels, layout rejection."""
import struct
import tempfile
import unittest
from pathlib import Path

import context

from mrm_quant import reader_adapter

SCHEMA = """<?xml version="1.0" encoding="utf-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema"
           xmlns:mstns="http://tempuri.org/XMLSchema.xsd"
           targetNamespace="http://tempuri.org/XMLSchema.xsd">
  <xs:element name="MSScanBin">
    <xs:complexType><xs:sequence>
      <xs:element minOccurs="0" maxOccurs="unbounded" name="ScanRecord" type="mstns:ScanRecordType"/>
    </xs:sequence></xs:complexType>
  </xs:element>
  <xs:complexType name="ScanRecordType"><xs:sequence>
%(fields)s
    <xs:element maxOccurs="unbounded" name="SpectrumParamValues" type="mstns:SpectrumParamsType"/>
  </xs:sequence></xs:complexType>
  <xs:complexType name="SpectrumParamsType"><xs:sequence>
    <xs:element name="SpectrumFormatID" type="xs:int"/>
    <xs:element name="SpectrumOffset" type="xs:long"/>
    <xs:element name="ByteCount" type="xs:int"/>
    <xs:element name="PointCount" type="xs:int"/>
    <xs:element name="MinY" type="xs:double"/>
    <xs:element name="MaxY" type="xs:double"/>
    <xs:element name="MinX" type="xs:double"/>
    <xs:element name="MaxX" type="xs:double"/>
  </xs:sequence></xs:complexType>
</xs:schema>
"""

# The scan record as declared by the datasets in hand: 3i d 2i 3d 4i 7d 2i.
RECORD_FIELDS = [('ScanID', 'int'), ('ScanMethodID', 'int'), ('TimeSegmentID', 'int'),
                 ('ScanTime', 'double'), ('MSLevel', 'int'), ('ScanType', 'int'),
                 ('TIC', 'double'), ('BasePeakMZ', 'double'), ('BasePeakValue', 'double'),
                 ('CycleNumber', 'int'), ('Status', 'int'), ('IonMode', 'int'),
                 ('IonPolarity', 'int'), ('Fragmentor', 'double'),
                 ('CollisionEnergy', 'double'), ('MzOfInterest', 'double'),
                 ('SamplingPeriod', 'double'), ('MeasuredMassRangeMin', 'double'),
                 ('MeasuredMassRangeMax', 'double'), ('Threshold', 'double'),
                 ('IsFragmentorDynamic', 'int'), ('IsCollisionEnergyDynamic', 'int')]

_RECORD = struct.Struct('<3i d 2i 3d 4i 7d 2i i q 2i 4d')
_DATA_OFFSET_POS = 0x58


def schema_text(extra_fields=()):
    fields = list(RECORD_FIELDS) + list(extra_fields)
    body = '\n'.join('    <xs:element name="%s" type="xs:%s"/>' % (name, kind)
                     for name, kind in fields)
    return SCHEMA % {'fields': body}


def synthetic_dataset(directory, *, byte_count=24, point_count=3, peak_bytes=24,
                      format_id=3, extra_fields=()):
    """A minimal .d tree with one MRM-like scan record and no method."""
    acqdata = Path(directory) / 'run.d' / 'AcqData'
    acqdata.mkdir(parents=True)
    (acqdata / 'MSScan.xsd').write_text(schema_text(extra_fields), encoding='utf-8')
    start = _DATA_OFFSET_POS + 4
    header = bytearray(start)
    struct.pack_into('<i', header, _DATA_OFFSET_POS, start)
    record = _RECORD.pack(1, 1, 1, 10.0, 2, 256, 6.0, 93.0, 3.0,
                          1, 0, 2, 0, 0.0, 30.0, 204.0, 0.0, 68.0, 93.0, 0.0,
                          0, 0, format_id, 0, byte_count, point_count,
                          0.0, 3.0, 68.0, 93.0)
    (acqdata / 'MSScan.bin').write_bytes(bytes(header) + record)
    (acqdata / 'MSPeak.bin').write_bytes(b'\x00' * peak_bytes)
    return acqdata.parent


class LayoutValidationTests(unittest.TestCase):
    def test_valid_synthetic_layout_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp)
            with reader_adapter.RunReader(path) as run:
                self.assertEqual(run.layout['record_stride_bytes'], 184)
                self.assertEqual(run.layout['schema_record_bytes'], 184)
                self.assertEqual(run.layout['scan_count'], 1)

    def test_byte_count_inconsistent_with_point_count_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp, byte_count=20, peak_bytes=20)
            with self.assertRaisesRegex(ValueError, 'unsupported_layout'):
                reader_adapter.RunReader(path)

    def test_peak_block_outside_file_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp, peak_bytes=8)
            with self.assertRaisesRegex(ValueError, 'unsupported_layout'):
                reader_adapter.RunReader(path)

    def test_trailing_unread_peak_bytes_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp, peak_bytes=64)
            with self.assertRaisesRegex(ValueError, 'unsupported_layout'):
                reader_adapter.RunReader(path)

    def test_unknown_spectrum_format_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp, format_id=7)
            with self.assertRaisesRegex(ValueError, 'unsupported_layout'):
                reader_adapter.RunReader(path)

    def test_schema_with_extra_field_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp, extra_fields=[('Reserved', 'double')])
            with self.assertRaisesRegex(ValueError, 'unsupported_layout'):
                reader_adapter.RunReader(path)

    def test_missing_schema_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = synthetic_dataset(tmp)
            (path / 'AcqData' / 'MSScan.xsd').unlink()
            with self.assertRaisesRegex(ValueError, 'unsupported_layout'):
                reader_adapter.RunReader(path)


class RealDatasetTests(unittest.TestCase):
    def run_reader(self, key):
        path = context.dataset_path(key)
        if not path.is_dir():
            self.skipTest('Local raw fixture unavailable: ' + key)
        run = reader_adapter.RunReader(path)
        self.addCleanup(run.close)
        return run

    def test_mrm_records_are_canonical(self):
        run = self.run_reader('mrm_primary')
        records = run.records(lambda head: head['ms_level'] == 2)
        self.assertTrue(records)
        record = records[0]
        self.assertEqual(record['scan_type'], 256)
        self.assertEqual(len(record['product_mz']), len(record['intensity']))
        self.assertEqual(record['product_mz'], sorted(record['product_mz']))
        self.assertEqual(record['product_mz'], record['declared_product_mz'])

    def test_mrm_channels_come_from_the_method(self):
        run = self.run_reader('mrm_primary')
        mrm = [c for c in run.channels if c['channel_type'] == 'mrm']
        self.assertTrue(mrm)
        self.assertTrue(all(c['frame_count'] > 0 for c in mrm))
        self.assertTrue(all(c['precursor_mz'] is not None for c in mrm))
        ms1 = [c for c in run.channels if c['channel_type'] == 'ms1']
        self.assertTrue(ms1)

    def test_full_scan_fixture_has_no_mrm_channel_or_frame(self):
        run = self.run_reader('ms1_reference')
        self.assertEqual([c for c in run.channels if c['channel_type'] == 'mrm'], [])
        self.assertEqual(run.metadata()['mrm_frame_count'], 0)
        self.assertEqual(run.records(lambda head: head['ms_level'] == 2), [])

    def test_different_methods_have_different_fingerprints(self):
        primary = self.run_reader('mrm_primary').method_fingerprint()
        secondary = self.run_reader('mrm_secondary').method_fingerprint()
        self.assertNotEqual(primary, secondary)


if __name__ == '__main__':
    unittest.main()

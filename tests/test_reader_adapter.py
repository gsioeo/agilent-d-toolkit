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
    def run_reader(self, relative):
        path = context.dataset_path(relative)
        if not path.is_dir():
            self.skipTest('Local raw dataset unavailable: ' + relative)
        run = reader_adapter.RunReader(path)
        self.addCleanup(run.close)
        return run

    def test_september_mrm_record_is_canonical(self):
        run = self.run_reader('TSJ-0907/LQ/STD_S1.d')
        records = run.records(lambda head: head['scan_id'] == 1)
        self.assertEqual(len(records), 1)
        record = records[0]
        self.assertEqual(record['scan_type'], 256)
        self.assertEqual(record['ms_level'], 2)
        self.assertEqual(record['scan_method_id'], 1)
        self.assertEqual(record['time_segment_id'], 1)
        self.assertEqual(record['polarity'], 0)
        self.assertEqual(record['precursor_mz'], 204.)
        self.assertEqual(record['collision_energy_ev'], 30.)
        self.assertAlmostEqual(record['rt_min'], 10.001733333333334)
        self.assertEqual(record['product_mz'], [68., 81., 93.])
        self.assertEqual(record['intensity'],
                         [0.11013031005859375, 5.686073303222656, 1121.11181640625])
        self.assertEqual(record['declared_product_mz'], [68., 81., 93.])

    def test_september_channels_come_from_the_method(self):
        run = self.run_reader('TSJ-0907/LQ/STD_S1.d')
        mrm = [c for c in run.channels if c['channel_type'] == 'mrm']
        self.assertEqual([c['product_mz'] for c in mrm], [93., 81., 68.])
        self.assertEqual({c['precursor_mz'] for c in mrm}, {204.})
        self.assertEqual({c['collision_energy_ev'] for c in mrm}, {30.})
        self.assertEqual({c['dwell_ms'] for c in mrm}, {100.})
        self.assertEqual({c['is_istd'] for c in mrm}, {False})
        self.assertEqual({c['compound_name'] for c in mrm}, {'Ses1', 'Ses2', 'Ses3'})
        self.assertEqual({c['frame_count'] for c in mrm}, {2454})
        ms1 = [c for c in run.channels if c['channel_type'] == 'ms1']
        self.assertEqual([(c['mz_low'], c['mz_high']) for c in ms1], [(50., 500.)])
        self.assertEqual(ms1[0]['frame_count'], 2454)

    def test_august_channels_differ_from_september(self):
        run = self.run_reader('20260819-tsj/TSJ-0819-1.d')
        products = {c['product_mz'] for c in run.channels if c['channel_type'] == 'mrm'}
        self.assertEqual(products, {189., 93., 69.})

    def test_july_has_no_mrm_channel_or_frame(self):
        run = self.run_reader('20260729/STD.d')
        self.assertEqual([c for c in run.channels if c['channel_type'] == 'mrm'], [])
        self.assertEqual(run.metadata()['mrm_frame_count'], 0)
        self.assertEqual(run.records(lambda head: head['ms_level'] == 2), [])

    def test_method_fingerprint_is_shared_within_a_batch_and_differs_between(self):
        september = [self.run_reader(f'TSJ-0907/LQ/STD_S{i}.d').method_fingerprint()
                     for i in (1, 2)]
        august = self.run_reader('20260819-tsj/TSJ-0819-1.d').method_fingerprint()
        self.assertEqual(len(set(september)), 1)
        self.assertNotIn(august, september)


if __name__ == '__main__':
    unittest.main()

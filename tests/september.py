"""Helpers for the September batch used by more than one test module.

The parameters below describe this one batch and were read off its own inspected
traces; they are test fixtures, not defaults for other data.
"""
import context

from mrm_quant import integration, reader_adapter, transitions

BATCH_ID = 'TSJ-0907-LQ'
ANALYTE_ID = 'Ses1'
QUANTIFIER_PRODUCT_MZ = 93.
EXPECTED_RT_MIN = 14.43
RT_TOLERANCE_MIN = 0.08
SEARCH_MARGIN_MIN = 0.25
FIND_OPTIONS = dict(min_relative_height=0.05, min_separation_min=0.1)
BASELINE = 'linear_endpoints'

TOP_CONCENTRATION = 1.0
DILUTION_STEP = 2.0
LEVELS = 10
CONCENTRATION_UNIT = 'mg/mL'

# The series spans 512-fold, so an unweighted fit is dominated by the top levels
# and its intercept (about 364) sits above the weakest sample areas, which would
# turn a real peak into a negative back-calculation. 1/x weighting keeps the low
# levels informative; it is a configured choice, not a selection made by
# comparing coefficients of determination.
WEIGHTING = '1/x'
INTERCEPT = 'free'


def concentrations():
    """S1 at the top concentration, halved for every following level."""
    return [TOP_CONCENTRATION / DILUTION_STEP ** index for index in range(LEVELS)]


def dataset(name):
    return context.dataset_path('TSJ-0907/LQ/%s.d' % name)


def trace(name, product_mz=QUANTIFIER_PRODUCT_MZ):
    path = dataset(name)
    with reader_adapter.RunReader(path) as run:
        channel = next(c for c in run.channels
                       if c['channel_type'] == 'mrm' and c['product_mz'] == product_mz)
        fingerprint = run.method_fingerprint()
    result = transitions.extract_run(path, transitions.selector_from_channel(channel))
    result['method_fingerprint'] = fingerprint
    return result


def chosen_peak(name, product_mz=QUANTIFIER_PRODUCT_MZ):
    extracted = trace(name, product_mz)
    peaks = integration.find_peaks(
        extracted['rt_min'], extracted['intensity'], baseline=BASELINE,
        rt_range=(EXPECTED_RT_MIN - SEARCH_MARGIN_MIN,
                  EXPECTED_RT_MIN + SEARCH_MARGIN_MIN), **FIND_OPTIONS)
    selection = integration.select_peak(peaks, expected_rt_min=EXPECTED_RT_MIN,
                                        tolerance_min=RT_TOLERANCE_MIN)
    selection['method_fingerprint'] = extracted['method_fingerprint']
    return selection


def calibration_rows():
    """One calibration row per standard, priced by the configured series."""
    rows = []
    for level, concentration in enumerate(concentrations(), start=1):
        name = 'STD_S%d' % level
        selection = chosen_peak(name)
        peak = selection['peak']
        rows.append({'run_id': name, 'batch_id': BATCH_ID, 'analyte_id': ANALYTE_ID,
                     'method_fingerprint': selection['method_fingerprint'],
                     'role': 'calibration', 'level_id': 'S%d' % level,
                     'concentration': concentration,
                     'concentration_unit': CONCENTRATION_UNIT,
                     'response': None if peak is None else peak['area'],
                     'include': True, 'exclusion_reason': None,
                     'peak_status': selection['status']})
    return rows

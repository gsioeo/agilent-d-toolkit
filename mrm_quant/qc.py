"""Qualifier ratio, blank check, result-status aggregation and output guards.

Nothing here invents a concentration. Each check returns a status; the
aggregation decides whether a back-calculated value may be reported at all, and
keeps the diagnostic numbers even when it may not. Absent checks are reported as
not evaluated rather than as passed.
"""
import math
from pathlib import Path

#: Qualifier ratio is defined as qualifier area over quantifier area, and the
#: tolerance as the relative deviation from the reference ratio.
RATIO_DEFINITION = 'r = A_qualifier / A_quantifier; deviation = abs(r / r_ref - 1)'

#: Statuses that forbid reporting a concentration.
BLOCKING_STATUSES = ('no_peak', 'ambiguous_peak', 'insufficient_points',
                     'non_positive_area', 'negative_backcalc',
                     'below_calibration_range', 'above_calibration_range',
                     'ion_ratio_fail', 'rt_mismatch', 'invalid_quantifier',
                     'blank_contamination', 'below_validated_loq',
                     'calibration_failed', 'quantification_failed',
                     'internal_standard_failed')


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) \
        and math.isfinite(float(value))


def check_qualifier(quant_area, qual_area, *, reference_ratio, relative_tolerance,
                    apex_delta_min, rt_tolerance_min):
    """Compare a qualifier transition with its quantifier.

    The reference ratio must come from qualified standards; a sample is never
    allowed to define its own reference.
    """
    if not _finite(reference_ratio) or float(reference_ratio) <= 0:
        raise ValueError('invalid_reference_ratio: %r' % reference_ratio)
    if not _finite(relative_tolerance) or float(relative_tolerance) < 0:
        raise ValueError('invalid_tolerance: %r' % relative_tolerance)
    result = {'status': 'ok', 'ratio': None, 'reference_ratio': float(reference_ratio),
              'relative_deviation': None, 'relative_tolerance': float(relative_tolerance),
              'apex_delta_min': apex_delta_min, 'rt_tolerance_min': rt_tolerance_min,
              'ratio_definition': RATIO_DEFINITION}
    if not _finite(quant_area) or float(quant_area) <= 0:
        result['status'] = 'invalid_quantifier'
        return result
    if not _finite(qual_area):
        result['status'] = 'invalid_qualifier'
        return result
    ratio = float(qual_area) / float(quant_area)
    result['ratio'] = ratio
    result['relative_deviation'] = abs(ratio / float(reference_ratio) - 1.)
    if apex_delta_min is not None and rt_tolerance_min is not None \
            and abs(float(apex_delta_min)) > float(rt_tolerance_min):
        result['status'] = 'rt_mismatch'
    elif result['relative_deviation'] > float(relative_tolerance):
        result['status'] = 'ion_ratio_fail'
    return result


def check_blank(blank_area, *, limit_area=None, blank_run_id=None):
    """Blank and carry-over check against a pre-set area limit."""
    result = {'status': 'not_evaluated', 'blank_area': blank_area,
              'limit_area': limit_area, 'blank_run_id': blank_run_id}
    if limit_area is None:
        return result
    if blank_area is None:
        result['status'] = 'no_peak'
        return result
    if not _finite(blank_area) or not _finite(limit_area):
        raise ValueError('nonfinite: blank area %r or limit %r' % (blank_area, limit_area))
    result['status'] = 'ok' if float(blank_area) <= float(limit_area) else 'blank_contamination'
    return result


def aggregate_status(*, peak_selection=None, integration=None, quantification=None,
                     qualifiers=(), blank=None, independent_qc=None,
                     loq_concentration=None, calibration_failure=None):
    """Final status and reportable concentration for one run and analyte.

    Returns ``status``, ``reported_concentration``, ``validated`` and the reasons
    behind them. A concentration is reported only when nothing blocks it; a
    missing quality-control check leaves ``validated`` false instead of claiming
    a pass.
    """
    reasons = []
    status = 'ok'

    def block(code, detail):
        nonlocal status
        if status == 'ok':
            status = code
        reasons.append(detail)

    if calibration_failure is not None:
        block('calibration_failed', 'calibration failed: %s' % calibration_failure)
    if peak_selection is not None and peak_selection.get('status') != 'ok':
        block(peak_selection['status'], 'peak selection: %s' % peak_selection['status'])
    if integration is not None and integration.get('status') not in (None, 'ok'):
        block(integration['status'], 'integration: %s' % integration['status'])
    if quantification is not None and quantification.get('status') != 'ok':
        detail = quantification.get('error') or quantification['status']
        block(quantification['status'], 'quantification: %s' % detail)
    for qualifier in qualifiers or ():
        if qualifier.get('status') != 'ok':
            block(qualifier['status'], 'qualifier %s: %s'
                  % (qualifier.get('channel_id', ''), qualifier['status']))
    if blank is not None and blank.get('status') == 'blank_contamination':
        block('blank_contamination', 'blank %s exceeds its limit'
              % (blank.get('blank_run_id') or ''))

    concentration = None
    if status == 'ok' and quantification is not None:
        concentration = quantification.get('original_concentration')
        if loq_concentration is not None and concentration is not None \
                and concentration < float(loq_concentration):
            status = 'below_validated_loq'
            reasons.append('below the validated limit of quantification')
            concentration = None

    validated = True
    if independent_qc is None or independent_qc.get('status') == 'not_evaluated':
        validated = False
        reasons.append('no independent quality control was evaluated')
    elif independent_qc.get('status') != 'ok':
        validated = False
        block('qc_fail', 'independent quality control: %s' % independent_qc['status'])
        concentration = None
    if blank is None or blank.get('status') == 'not_evaluated':
        validated = False
        reasons.append('no blank limit was evaluated')

    return {'status': status, 'reported_concentration': concentration,
            'validated': validated, 'reasons': reasons,
            'vial_concentration': None if quantification is None
            else quantification.get('vial_concentration')}


def validate_output_path(path, *, raw_datasets=()):
    """Refuse an output directory that could damage raw data or earlier results.

    Rejects a path inside any raw dataset, any path below a ``.d`` directory, and
    any path that already exists, after resolving symbolic links.
    """
    candidate = Path(path)
    resolved = candidate.expanduser().resolve()
    for dataset in raw_datasets or ():
        dataset_resolved = Path(dataset).expanduser().resolve()
        if resolved == dataset_resolved or dataset_resolved in resolved.parents:
            raise ValueError('unsafe_output: %s is inside the raw dataset %s'
                             % (resolved, dataset_resolved))
    for part in (resolved,) + tuple(resolved.parents):
        if part.name.lower().endswith('.d'):
            raise ValueError('unsafe_output: %s is inside a .d directory' % resolved)
    if candidate.exists() or resolved.exists():
        raise ValueError('unsafe_output: %s already exists; results are never '
                         'written into an existing directory' % resolved)
    return resolved

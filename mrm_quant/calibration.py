"""Same-batch calibration fitting and back-calculation.

The model is the one the method declares in advance: a straight line
``A = a * C_vial + b`` fitted by weighted least squares, with the weighting and
the intercept chosen by configuration rather than by whichever variant produces
the best coefficient of determination. Only calibration standards of the same
batch, analyte, method and concentration unit take part. Blanks, quality
controls and unknowns never enter the fit, replicates count as points but not as
levels, and no point is dropped without a recorded reason.
"""
import math

WEIGHTINGS = ('none', '1/x', '1/x2')
INTERCEPT_MODES = ('free', 'zero')
INVERSE_WEIGHTINGS = ('1/x', '1/x2')
RESPONSE_MODES = ('external', 'internal')

#: Fewest distinct concentrations this implementation will fit.
MIN_LEVELS = 3

WEIGHT_DEFINITIONS = {'none': 'w = 1', '1/x': 'w = 1 / C', '1/x2': 'w = 1 / C^2'}


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) \
        and math.isfinite(float(value))


def _weight(concentration, weighting):
    if weighting == 'none':
        return 1.
    if concentration <= 0:
        raise ValueError('zero_weight_domain: weighting %s needs a positive '
                         'concentration, got %r' % (weighting, concentration))
    if weighting == '1/x':
        return 1. / concentration
    return 1. / (concentration * concentration)


def _single(values, code, label):
    unique = {value for value in values}
    if len(unique) > 1:
        raise ValueError('%s: calibration points span %s %s'
                         % (code, label, sorted(map(repr, unique))))
    return unique.pop() if unique else None


def _fit_points(rows, weighting):
    points = []
    excluded = []
    for row in rows:
        if (row.get('role') or '').strip().lower() != 'calibration':
            continue
        if not row.get('include', True):
            reason = row.get('exclusion_reason')
            if not reason:
                raise ValueError('exclusion_reason_required: run %r is excluded '
                                 'without a reason' % row.get('run_id'))
            excluded.append({'run_id': row.get('run_id'), 'reason': reason})
            continue
        concentration = row.get('concentration')
        response = row.get('response')
        if concentration is None:
            raise ValueError('missing_concentration: run %r has no concentration'
                             % row.get('run_id'))
        if not _finite(concentration):
            raise ValueError('nonfinite: run %r concentration is %r'
                             % (row.get('run_id'), concentration))
        if concentration < 0:
            raise ValueError('negative_concentration: run %r concentration is %r'
                             % (row.get('run_id'), concentration))
        if response is None:
            raise ValueError('missing_response: run %r has no response'
                             % row.get('run_id'))
        if not _finite(response):
            raise ValueError('nonfinite: run %r response is %r'
                             % (row.get('run_id'), response))
        points.append({'run_id': row.get('run_id'), 'level_id': row.get('level_id'),
                       'concentration': float(concentration),
                       'response': float(response),
                       'weight': _weight(float(concentration), weighting)})
    return points, excluded


def fit_calibration(rows, *, weighting, intercept):
    """Weighted linear calibration over the calibration standards in *rows*."""
    if weighting not in WEIGHTINGS:
        raise ValueError('invalid_weighting: %r is not one of %s'
                         % (weighting, ', '.join(WEIGHTINGS)))
    if intercept not in INTERCEPT_MODES:
        raise ValueError('invalid_intercept: %r is not one of %s'
                         % (intercept, ', '.join(INTERCEPT_MODES)))

    calibration_rows = [row for row in rows
                        if (row.get('role') or '').strip().lower() == 'calibration']
    batch_id = _single([row.get('batch_id') for row in calibration_rows],
                       'batch_mismatch', 'batch')
    analyte_id = _single([row.get('analyte_id') for row in calibration_rows],
                         'analyte_mismatch', 'analyte')
    fingerprint = _single([row.get('method_fingerprint') for row in calibration_rows],
                          'method_mismatch', 'method fingerprint')
    unit = _single([row.get('concentration_unit') for row in calibration_rows],
                   'unit_mismatch', 'concentration unit')

    points, excluded = _fit_points(rows, weighting)
    if len(points) < 2:
        raise ValueError('insufficient_points: %d usable calibration point(s)'
                         % len(points))
    levels = sorted({point['concentration'] for point in points})
    if len(levels) < MIN_LEVELS:
        raise ValueError('insufficient_levels: %d distinct concentration(s), '
                         'at least %d are required' % (len(levels), MIN_LEVELS))

    sum_w = sum(p['weight'] for p in points)
    sum_wx = sum(p['weight'] * p['concentration'] for p in points)
    sum_wy = sum(p['weight'] * p['response'] for p in points)
    sum_wxx = sum(p['weight'] * p['concentration'] ** 2 for p in points)
    sum_wxy = sum(p['weight'] * p['concentration'] * p['response'] for p in points)

    if intercept == 'zero':
        if sum_wxx <= 0:
            raise ValueError('degenerate_fit: all concentrations are zero')
        slope = sum_wxy / sum_wxx
        offset = 0.
    else:
        denominator = sum_w * sum_wxx - sum_wx * sum_wx
        if denominator == 0 or not math.isfinite(denominator):
            raise ValueError('degenerate_fit: the design matrix is singular')
        slope = (sum_w * sum_wxy - sum_wx * sum_wy) / denominator
        offset = (sum_wy - slope * sum_wx) / sum_w
    if not math.isfinite(slope) or not math.isfinite(offset):
        raise ValueError('nonfinite: the fit produced %r and %r' % (slope, offset))
    if slope <= 0:
        raise ValueError('non_positive_slope: slope is %r' % slope)

    residuals, bias, back = [], [], []
    for point in points:
        predicted = slope * point['concentration'] + offset
        residuals.append(point['response'] - predicted)
        recovered = (point['response'] - offset) / slope
        back.append(recovered)
        bias.append(None if point['concentration'] == 0
                    else (recovered - point['concentration']) / point['concentration'] * 100.)

    mean_response = sum(p['response'] for p in points) / len(points)
    total = sum((p['response'] - mean_response) ** 2 for p in points)
    residual = sum(value * value for value in residuals)
    r2 = None if total == 0 else 1. - residual / total

    weighted_mean = sum(p['weight'] * p['response'] for p in points) / sum_w
    weighted_total = sum(p['weight'] * (p['response'] - weighted_mean) ** 2 for p in points)
    weighted_residual = sum(point['weight'] * value * value
                            for point, value in zip(points, residuals))
    r2_weighted = None if weighted_total == 0 else 1. - weighted_residual / weighted_total

    for point, recovered, deviation, resid in zip(points, back, bias, residuals):
        point['back_calculated_concentration'] = recovered
        point['bias_pct'] = deviation
        point['residual'] = resid

    equation = 'A = %.6g * C %s %.6g' % (slope, '-' if offset < 0 else '+', abs(offset))

    return {'slope': slope, 'intercept': offset, 'weighting': weighting,
            'weight_definition': WEIGHT_DEFINITIONS[weighting],
            'intercept_mode': intercept, 'range': [levels[0], levels[-1]],
            'n_levels': len(levels), 'n_points': len(points),
            'equation': equation, 'r2': r2, 'r2_weighted': r2_weighted,
            'r2_definition': ('r2 = 1 - sum((A - fit)^2) / sum((A - mean(A))^2); '
                              'r2_weighted applies the fit weights to both sums. '
                              'Descriptive only: values from different weightings '
                              'are not comparable and do not select a model.'),
            'backcalc_bias_pct': bias, 'residuals': residuals,
            'r2_unweighted': r2, 'excluded_run_ids': [item['run_id'] for item in excluded],
            'exclusions': excluded, 'batch_id': batch_id, 'analyte_id': analyte_id,
            'method_fingerprint': fingerprint, 'concentration_unit': unit,
            'points': points, 'model': 'linear'}


def response_value(area, *, mode, is_area=None):
    """Quantification response: the area itself, or the ratio to a fixed internal standard."""
    if mode not in RESPONSE_MODES:
        raise ValueError('invalid_response_mode: %r is not one of %s'
                         % (mode, ', '.join(RESPONSE_MODES)))
    if area is None or not _finite(area):
        raise ValueError('invalid_area: %r' % area)
    if mode == 'external':
        return float(area)
    if is_area is None or not _finite(is_area) or float(is_area) <= 0:
        raise ValueError('invalid_internal_standard: %r is not a positive area' % is_area)
    return float(area) / float(is_area)


def quantify(response, model, *, dilution_factor, batch_id, method_fingerprint):
    """Back-calculate one concentration from a response and a calibration model.

    Only the arithmetic and the calibration range are decided here. The reported
    result still depends on the quality-control aggregation.
    """
    if model.get('batch_id') != batch_id:
        raise ValueError('batch_mismatch: model belongs to batch %r, not %r'
                         % (model.get('batch_id'), batch_id))
    if model.get('method_fingerprint') != method_fingerprint:
        raise ValueError('method_mismatch: model belongs to method %r, not %r'
                         % (model.get('method_fingerprint'), method_fingerprint))
    if dilution_factor is None or not _finite(dilution_factor) or float(dilution_factor) <= 0:
        raise ValueError('invalid_dilution: %r is not a positive factor' % dilution_factor)
    if response is None or not _finite(response):
        raise ValueError('nonfinite: response is %r' % response)

    slope = model['slope']
    offset = model['intercept']
    vial = (float(response) - offset) / slope
    low, high = model['range']
    result = {'response': float(response), 'vial_concentration': vial,
              'original_concentration': None, 'dilution_factor': float(dilution_factor),
              'concentration_unit': model.get('concentration_unit'),
              'batch_id': batch_id, 'method_fingerprint': method_fingerprint,
              'calibration_range': [low, high], 'status': 'ok'}
    if vial < 0:
        result['status'] = 'negative_backcalc'
    elif vial < low:
        result['status'] = 'below_calibration_range'
    elif vial > high:
        result['status'] = 'above_calibration_range'
    else:
        result['original_concentration'] = vial * float(dilution_factor)
    return result

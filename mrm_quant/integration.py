"""Peak integration inside an expected retention-time window.

Areas are plain trapezoids over the acquired, unevenly spaced samples, in
minutes, with either no baseline or a straight line between the two window
endpoints. Window endpoints that fall between samples are interpolated; nothing
is extrapolated beyond the acquired data. Residuals after baseline subtraction
keep their sign, so a sloping baseline does not bias the area upwards.
"""
import math

AREA_UNIT = 'stored_intensity*min'
BASELINE_MODES = ('none', 'linear_endpoints')

#: Comparison slack for window edges landing exactly on a sample time.
EDGE_EPS = 1e-12


def _check_baseline(baseline):
    if baseline not in BASELINE_MODES:
        raise ValueError('invalid_baseline: %r is not one of %s'
                         % (baseline, ', '.join(BASELINE_MODES)))


def _validate_slice(rt, signal, start, stop, max_gap_min):
    """Validate the samples this integration will consult, and nothing else."""
    previous = None
    for index in range(start, stop):
        time = rt[index]
        value = signal[index]
        if time is None or not math.isfinite(time):
            raise ValueError('nonfinite: sample %d has no usable time' % index)
        if value is None:
            raise ValueError('missing_data: sample %d has no intensity' % index)
        value = float(value)
        if not math.isfinite(value):
            raise ValueError('nonfinite: sample %d intensity is %r' % (index, value))
        if previous is not None:
            if time <= previous:
                raise ValueError('non_increasing_time: sample %d repeats or precedes %r'
                                 % (index, previous))
            if max_gap_min is not None and (time - previous) > max_gap_min:
                raise ValueError('acquisition_gap: %g min between samples %d and %d '
                                 'exceeds %g min' % (time - previous, index - 1,
                                                     index, max_gap_min))
        previous = time


def _interpolate(rt, signal, target, lower, upper):
    span = rt[upper] - rt[lower]
    if span <= 0:
        raise ValueError('non_increasing_time: samples %d and %d share a time'
                         % (lower, upper))
    weight = (target - rt[lower]) / span
    return float(signal[lower]) + weight * (float(signal[upper]) - float(signal[lower]))


def _window_points(rt, signal, bounds_min, max_gap_min):
    low, high = float(bounds_min[0]), float(bounds_min[1])
    if not (math.isfinite(low) and math.isfinite(high)):
        raise ValueError('nonfinite: integration bounds are %r' % (bounds_min,))
    if high <= low:
        raise ValueError('invalid_bounds: %g is not below %g' % (low, high))
    if len(rt) != len(signal):
        raise ValueError('length_mismatch: %d times for %d intensities'
                         % (len(rt), len(signal)))
    if len(rt) < 2:
        raise ValueError('insufficient_data: at least two samples are required')
    if low < rt[0] - EDGE_EPS or high > rt[-1] + EDGE_EPS:
        raise ValueError('window_outside_data: [%g, %g] is not covered by [%r, %r]'
                         % (low, high, rt[0], rt[-1]))

    inside = [i for i, time in enumerate(rt)
              if time is not None and low - EDGE_EPS <= time <= high + EDGE_EPS]
    if not inside:
        bracket = [i for i, time in enumerate(rt) if time is not None and time < low]
        start = bracket[-1] if bracket else 0
        stop = min(start + 2, len(rt))
    else:
        start = inside[0] - 1 if inside[0] > 0 and rt[inside[0]] > low + EDGE_EPS else inside[0]
        stop = inside[-1] + 2 if inside[-1] + 1 < len(rt) and rt[inside[-1]] < high - EDGE_EPS \
            else inside[-1] + 1
    _validate_slice(rt, signal, start, stop, max_gap_min)

    points = []
    if inside and abs(rt[inside[0]] - low) <= EDGE_EPS:
        pass
    else:
        upper = inside[0] if inside else stop - 1
        points.append((low, _interpolate(rt, signal, low, upper - 1, upper)))
    for index in inside:
        points.append((rt[index], float(signal[index])))
    if inside and abs(rt[inside[-1]] - high) <= EDGE_EPS:
        pass
    else:
        lower = inside[-1] if inside else start
        points.append((high, _interpolate(rt, signal, high, lower, lower + 1)))
    points.sort(key=lambda item: item[0])
    return points


def integrate_peak(rt, signal, *, bounds_min, baseline='linear_endpoints',
                   max_gap_min=None, min_points=3):
    """Trapezoidal area between ``bounds_min``, in ``stored_intensity*min``.

    ``baseline='none'`` integrates the signal as stored; ``'linear_endpoints'``
    subtracts the straight line joining the two window endpoints. A window with
    fewer than ``min_points`` samples, or a non-positive net area, returns a
    status other than ``ok`` and is not a valid quantification peak.
    """
    _check_baseline(baseline)
    points = _window_points(rt, signal, bounds_min, max_gap_min)
    result = {'status': 'ok', 'area': None, 'area_unit': AREA_UNIT,
              'baseline': baseline, 'bounds_min': [float(bounds_min[0]), float(bounds_min[1])],
              'point_count': len(points), 'apex_rt_min': None, 'apex_intensity': None,
              'height_above_baseline': None}
    if len(points) < min_points:
        result['status'] = 'insufficient_points'
        return result

    if baseline == 'linear_endpoints':
        (t0, y0), (t1, y1) = points[0], points[-1]
        slope = (y1 - y0) / (t1 - t0)

        def base(time):
            return y0 + slope * (time - t0)
    else:
        def base(time):
            return 0.

    area = 0.
    apex_index = 0
    apex_excess = None
    for index, (time, value) in enumerate(points):
        excess = value - base(time)
        if apex_excess is None or excess > apex_excess:
            apex_excess, apex_index = excess, index
        if index:
            previous_time, previous_value = points[index - 1]
            previous_excess = previous_value - base(previous_time)
            area += 0.5 * (excess + previous_excess) * (time - previous_time)

    result['area'] = area
    result['apex_rt_min'] = points[apex_index][0]
    result['apex_intensity'] = points[apex_index][1]
    result['height_above_baseline'] = apex_excess
    if not math.isfinite(area):
        raise ValueError('nonfinite: integrated area is %r' % area)
    if area <= 0:
        result['status'] = 'non_positive_area'
    return result


def find_peaks(rt, signal, *, baseline='linear_endpoints', max_gap_min=None,
               rt_range=None, min_relative_height=0.02, min_absolute_height=None,
               boundary_fraction=0.05, min_separation_min=0.04, max_peaks=None,
               min_points=3):
    """Candidate peaks with apex and area, for :func:`select_peak` to choose from.

    Local maxima are traced down to the neighbouring minima and cut at
    ``boundary_fraction`` of the apex height, then integrated with the same
    rules as :func:`integrate_peak`. Nothing here knows the analyte; the
    retention-time decision belongs to :func:`select_peak`.
    """
    usable = [i for i, (time, value) in enumerate(zip(rt, signal))
              if time is not None and value is not None
              and math.isfinite(time) and math.isfinite(float(value))
              and (rt_range is None or rt_range[0] <= time <= rt_range[1])]
    if len(usable) < 3:
        return []
    heights = [float(signal[i]) for i in usable]
    ceiling = max(heights)
    if ceiling <= 0:
        return []
    floor = min_absolute_height if min_absolute_height is not None \
        else min_relative_height * ceiling

    candidates = []
    for position in range(1, len(usable) - 1):
        index = usable[position]
        value = float(signal[index])
        if value < floor:
            continue
        if value > float(signal[usable[position - 1]]) and value >= float(signal[usable[position + 1]]):
            candidates.append(position)
    candidates.sort(key=lambda position: float(signal[usable[position]]), reverse=True)

    peaks, taken = [], []
    for position in candidates:
        if max_peaks is not None and len(peaks) >= max_peaks:
            break
        apex = usable[position]
        if any(abs(rt[apex] - rt[other]) < min_separation_min for other in taken):
            continue
        cut = float(signal[apex]) * boundary_fraction
        left = position
        while left > 0 and float(signal[usable[left - 1]]) < float(signal[usable[left]]) \
                and float(signal[usable[left]]) > cut:
            left -= 1
        right = position
        while right < len(usable) - 1 \
                and float(signal[usable[right + 1]]) < float(signal[usable[right]]) \
                and float(signal[usable[right]]) > cut:
            right += 1
        bounds = (rt[usable[left]], rt[usable[right]])
        if bounds[1] <= bounds[0]:
            continue
        try:
            integrated = integrate_peak(rt, signal, bounds_min=bounds, baseline=baseline,
                                        max_gap_min=max_gap_min, min_points=min_points)
        except ValueError as error:
            peaks.append({'apex_rt_min': rt[apex], 'apex_intensity': float(signal[apex]),
                          'area': None, 'area_unit': AREA_UNIT, 'bounds_min': list(bounds),
                          'baseline': baseline, 'status': str(error).split(':', 1)[0]})
            taken.append(apex)
            continue
        peak = {'apex_rt_min': rt[apex], 'apex_intensity': float(signal[apex]),
                'area': integrated['area'], 'area_unit': AREA_UNIT,
                'bounds_min': list(bounds), 'baseline': baseline,
                'height_above_baseline': integrated['height_above_baseline'],
                'status': integrated['status'], 'point_count': integrated['point_count']}
        peaks.append(peak)
        taken.append(apex)
    return sorted(peaks, key=lambda peak: peak['apex_rt_min'])


def select_peak(peaks, *, expected_rt_min, tolerance_min):
    """The single candidate inside the expected retention-time window.

    No candidate is ``no_peak``; more than one is ``ambiguous_peak`` and returns
    no peak for review. The most intense peak of the run is never chosen on
    intensity alone.
    """
    if tolerance_min is None or tolerance_min < 0:
        raise ValueError('invalid_tolerance: %r' % tolerance_min)
    inside = [peak for peak in peaks
              if peak.get('apex_rt_min') is not None
              and abs(float(peak['apex_rt_min']) - float(expected_rt_min)) <= tolerance_min]
    result = {'status': 'ok', 'peak': None, 'candidates': len(inside),
              'expected_rt_min': float(expected_rt_min),
              'tolerance_min': float(tolerance_min)}
    if not inside:
        result['status'] = 'no_peak'
    elif len(inside) > 1:
        result['status'] = 'ambiguous_peak'
    else:
        result['peak'] = inside[0]
    return result

"""Peak integration inside an expected retention-time window."""


def integrate_peak(rt, signal, *, bounds_min, baseline, max_gap_min):
    raise NotImplementedError('integration.integrate_peak is implemented in task 4')


def find_peaks(rt, signal, **options):
    raise NotImplementedError('integration.find_peaks is implemented in task 4')


def select_peak(peaks, *, expected_rt_min, tolerance_min):
    raise NotImplementedError('integration.select_peak is implemented in task 4')

"""Same-batch calibration fitting and back-calculation."""


def fit_calibration(rows, *, weighting, intercept):
    raise NotImplementedError('calibration.fit_calibration is implemented in task 5')


def response_value(area, *, mode, is_area=None):
    raise NotImplementedError('calibration.response_value is implemented in task 5')


def quantify(response, model, *, dilution_factor, batch_id, method_fingerprint):
    raise NotImplementedError('calibration.quantify is implemented in task 5')

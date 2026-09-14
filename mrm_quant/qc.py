"""Qualifier ratio, blank and result-status aggregation, plus output guards."""


def check_qualifier(quant_area, qual_area, *, reference_ratio, relative_tolerance,
                    apex_delta_min, rt_tolerance_min):
    raise NotImplementedError('qc.check_qualifier is implemented in task 6')


def check_blank(blank_area, *, limit_area):
    raise NotImplementedError('qc.check_blank is implemented in task 6')


def aggregate_status(**parts):
    raise NotImplementedError('qc.aggregate_status is implemented in task 6')


def validate_output_path(path, *, raw_datasets):
    raise NotImplementedError('qc.validate_output_path is implemented in task 6')

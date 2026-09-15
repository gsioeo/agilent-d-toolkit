"""mrm_quant: MRM transition extraction and same-batch calibration quantification.

Read-only with respect to the raw Agilent ``.d`` datasets and to the existing
``ingest/`` scripts. Every analytical parameter comes from configuration files;
nothing about a particular batch, transition or concentration is hardcoded here.
"""

VERSION = 'v2'

__all__ = ['VERSION']

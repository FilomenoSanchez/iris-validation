"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

from iris_validation._defs import *
from iris_validation import metrics, interface


def generate_report(latest_model_path,
                    previous_model_path=None,
                    latest_reflections_path=None,
                    previous_reflections_path=None,
                    output_dir='./iris_report_output/',
                    mode='',
                    sequence_path = None,
                    distance_prediction_format=None,
                    distance_prediction_path=None,
                    map_align_exe='map_align',
                    dssp_exe='mkdssp'):

    latest_metrics_model = metrics.generate_metrics_model(latest_model_path, latest_reflections_path,
                                                          f_seq=sequence_path, f_distpred=distance_prediction_path,
                                                          distpred_format=distance_prediction_format,
                                                          map_align_exe=map_align_exe, dssp_exe=dssp_exe)
    previous_metrics_model = None
    if previous_model_path is not None:
        previous_metrics_model = metrics.generate_metrics_model(previous_model_path, previous_reflections_path,
                                                                f_seq=sequence_path, f_distpred=distance_prediction_path,
                                                                distpred_format=distance_prediction_format,
                                                                map_align_exe=map_align_exe, dssp_exe=dssp_exe)

    interface.build_report(latest_metrics_model, previous_metrics_model, output_dir, mode)

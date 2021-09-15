"""
Copyright 2020 William Rochira at York Structural Biology Laboratory

- Generates an example Iris reoprt for a given PDB code using model and reflections
  data from PDB-REDO
"""

import os
import sys
import time

from iris_validation import generate_report
from iris_tools.common import get_from_pdb_redo


INPUT_DIR = 'example_input'
OUTPUT_DIR_PREFIX = 'example_report'

if __name__ == '__main__':
    #pdb_id = raw_input('Enter PDB code: ')

    latest_model_path = '/home/filo/tmp/latest_pdb.pdb'
    previous_model_path = '/home/filo/tmp/previous_pdb.pdb'
    latest_reflections_path = '/home/filo/Documents/Register_errors/5MSZ/5msz_phases.mtz'
    previous_reflections_path = '/home/filo/Documents/Register_errors/5MSZ/5msz_phases.mtz'
    sequence_path = '/home/filo/Documents/Register_errors/5MSZ/5msz.fasta'
    distpred_path = '/home/filo/Documents/Register_errors/5MSZ/TR049386_results/seq.npz'
    distpred_format = 'rosettanpz'
    map_align_exe = '/home/filo/opt/map_align/map_align'

    all_files_present = True
    for path in (latest_model_path, latest_reflections_path):
        if not os.path.exists(path):
            all_files_present = False
            continue

    if not all_files_present:
       print('Failed to get files.')
       exit(1)

    generate_report(latest_model_path,
                    previous_model_path,
                    latest_reflections_path,
                    previous_reflections_path,
                    sequence_path=sequence_path,
                    distance_prediction_path=distpred_path,
                    distance_prediction_format=distpred_format,
                    output_dir='/home/filo/iris_output',
                    mode='', map_align_exe=map_align_exe)

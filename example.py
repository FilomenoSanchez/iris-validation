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
OUTPUT_DIR_PREFIX = 'example_report_{}'

if __name__ == '__main__':

    pdb_id = raw_input('Enter PDB code: ')

    latest_model_path = os.path.join(INPUT_DIR, '{}_final.pdb'.format(pdb_id))
    previous_model_path = os.path.join(INPUT_DIR, '{}_0cyc.pdb'.format(pdb_id))
    latest_reflections_path = os.path.join(INPUT_DIR, '{}_final.mtz'.format(pdb_id))
    previous_reflections_path = os.path.join(INPUT_DIR, '{}_0cyc.mtz'.format(pdb_id))
    sequence_path = os.path.join(INPUT_DIR, '{}.fasta'.format(pdb_id))
    distpred_path = os.path.join(INPUT_DIR, '{}.npz'.format(pdb_id))
    distpred_format = 'rosettanpz'
    map_align_exe = 'map_align'
    dssp_exe = 'mkdssp'

    all_files_present = True
    for path in (latest_model_path, latest_reflections_path, sequence_path, distpred_path):
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
                    output_dir=OUTPUT_DIR_PREFIX.format(pdb_id),
                    mode='', map_align_exe=map_align_exe, dssp_exe=dssp_exe)

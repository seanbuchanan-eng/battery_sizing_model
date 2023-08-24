"""
This file loads data one channel at a time and saves the created objects to disk.

Processing the data in this way ensures that there is enough RAM for the computations.
"""

import os
os.chdir('..')
print(os.getcwd())
import pickle
import numpy as np
from copy import deepcopy
import primed_data_processing.primed_utils as primed_utils
from primed_data_processing.cellbuilder import CellBuilder
from primed_data_processing.arbin_cycler import ArbinBatch, ArbinCell

def load_B6_test(cell_builder: CellBuilder, 
               prepath: str,
               folder1: str,
               folder2: str,
               channel_numbers: list | tuple, 
               cell_numbers: list | tuple, 
               steps: dict[str: list[int]],
               insert_step: int
               ) -> ArbinBatch:

    # list for holding processed cells
    arbin_cells = []

    # make new steps for T15
    first_test_steps = steps
    second_test_steps = deepcopy(steps)
    second_test_steps['degradation'] = list(np.array(second_test_steps['degradation']) + 1)
    second_test_steps['degradation'].insert(0,insert_step)

    # loop over channel numbers
    for channel_idx, channel in enumerate(channel_numbers):
        print(f'Processing channel {channel}')

        # append new cell to cells processed cells list
        arbin_cells.append(ArbinCell(cell_numbers[channel_idx], channel))

        # make subfolder name in raws folder
        # two folders are needed because the tests were modified part way through
        # hence two seperate folders contain the "different" data.
        first_test_folder_name = folder1 + f'/Channel_{channel}/'
        second_test_folder_name = folder2 + f'/Channel_{channel}/'

        # get directory of the current folder
        first_test_directory = os.fsencode(prepath+first_test_folder_name)
        second_test_directory = os.fsencode(prepath+second_test_folder_name)

        # sort directory into chronological order for read_B6_csv_data()
        first_test_sorted_dir = sorted(os.listdir(first_test_directory), key=lambda file: int(os.fsdecode(file).split('.')[1]))
        second_test_sorted_dir = sorted(os.listdir(second_test_directory), key=lambda file: int(os.fsdecode(file).split('.')[1]))

        # load the files using cellbuilder
        primed_utils.load_files_from_dir(cell_builder, arbin_cells, first_test_sorted_dir, prepath, first_test_folder_name, first_test_steps, channel_idx)
        primed_utils.load_files_from_dir(cell_builder, arbin_cells, second_test_sorted_dir, prepath, second_test_folder_name, second_test_steps, channel_idx)

    # load cells into a batch
    return ArbinBatch(cells=arbin_cells)


def read_data(channel_numbers: list | tuple, 
              cell_numbers: list | tuple, 
              cell_builder: CellBuilder, 
              raws_prepath: str, 
              test_number: int
              ) -> None:
    for idx, channel in enumerate(channel_numbers):
        if test_number == 10:
            folder1 = "B6T10V0_1_2_3_4_9_10_11_12_13_14_15_16"
            folder2 = "B6T15V0_1_2_3_4_9_10_11_12_13_14_15_16"
            batch = load_B6_test(cell_builder, raws_prepath, folder1, folder2, channel_numbers[idx:idx+1], cell_numbers[idx:idx+1], steps, 21)
            with open(f"data_objects/B6T1015C{cell_numbers[idx]}_active_steps.pkl", "wb") as f:
                pickle.dump(batch[:], f)
        elif test_number == 11:
            folder1 = "B6T11V0_6_7"
            folder2 = "B6T16V0_6_7"
            batch = load_B6_test(cell_builder, raws_prepath, folder1, folder2, channel_numbers[idx:idx+1], cell_numbers[idx:idx+1], steps, 19)
            with open(f"data_objects/B6T1116C{cell_numbers[idx]}_active_steps.pkl", "wb") as f:
                pickle.dump(batch[:], f)
        elif test_number == 12:
            folder1 = "B6T12V0_5_8"
            folder2 = "B6T17V0_5_8"
            batch = load_B6_test(cell_builder, raws_prepath, folder1, folder2, channel_numbers[idx:idx+1], cell_numbers[idx:idx+1], steps, 19)
            with open(f"data_objects/B6T1217C{cell_numbers[idx]}_active_steps.pkl", "wb") as f:
                pickle.dump(batch[:], f)
        else:
            print("invalid test number.")
            return
        del batch

if __name__ == '__main__':
    
    cell_builder = CellBuilder()
    raws_prepath = 'C:/Users/seanb/OneDrive/Documents/PRIMED/export/batdat/MTC/raws/'

    # all channel and cell numbers from B6T10 and T15 in order
    channel_numbers = (1,2,3,4,9,10,11,12,13,14,15,16)
    cell_numbers = (9,10,11,12,1,2,3,4,5,6,7,8)
    # test 10 steps
    steps = {'characterization': [1,3,4,6,7,8,10,11,12,13,15,17,19,20],
            'degradation': [22,23,25,26,27,28,29,31]}

    read_data(channel_numbers, cell_numbers, cell_builder, raws_prepath, 10)

    # test 11
    channel_numbers = (6,7)
    cell_numbers = (14,15)
    steps = {'characterization': [1,3,4,6,7,8,10,11,12,13,15,17,18],
            'degradation': [20,21,22,23,25]}
    
    read_data(channel_numbers, cell_numbers, cell_builder, raws_prepath, 11)

    #test 12
    channel_numbers = (5,8)
    cell_numbers = (13,16)
    steps = {'characterization': [1,3,4,6,7,8,10,11,12,13,15,17,18],
            'degradation': [20]}
    
    read_data(channel_numbers, cell_numbers, cell_builder, raws_prepath, 12)
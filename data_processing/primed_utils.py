"""
Author: Sean Buchanan

This module contains useful functions for processing Arbin and Gamry battery data.

These functions come up continuously in practice from users so they have been 
added to a useful common location.
"""

import os
from primed_data_processing.cellbuilder import CellBuilder
from primed_data_processing.arbin_cycler import ArbinBatch, ArbinCell, ArbinStep
from primed_data_processing.gamry_eis import EisSweep, EisCycle, EisCell

#################################### ARBIN ######################################################

def load_files_from_dir(cell_builder: CellBuilder, 
                        arbin_cells: list[ArbinCell], 
                        sorted_dir: list[bytes], 
                        prepath: str, 
                        folder_name: str, 
                        steps: dict[str: list[int]], 
                        channel_idx: int
                        ) -> None:
    """
    Load all csv files in a directory.

    Parameters
    ----------
    ``cell_builder`` \: ``CellBuilder``
        ``CellBuilder`` used to access data reading methods.
    ``arbin_cells`` \: ``list[ArbinCell]``
        list of arbin cells to add data to.
    ``sorted_dir`` \: ``list[bytes]``
        list of files in the directory sorted by chronological order.
    ``prepath`` \: ``str``
        path to the folder containing all battery data.
    ``folder_name`` \: ``str``
        name of the folder to be read.
    ``steps`` \: ``dict[str: list[int]]``
        dictionary of steps to be read.
    ``channel_idx`` \: ``int``
        Which channel to load data from.
    """

    # loop over all files in the current directory
    for file in sorted_dir:
        # get filename
        filename = os.fsdecode(file)

        # ignore .xlsx files
        if filename.endswith('.csv'):
            # process file with CellBuilder method
            cell_builder.read_B6_csv_data(
                arbin_cells[channel_idx], # current cell being processed
                prepath+folder_name+filename, # path to file being processed
                steps, # get OCV charge and discharg steps
                verbose=False # minimal printouts
            )

def load_B6T10(cell_builder: CellBuilder, 
               prepath: str, 
               channel_numbers: list | tuple, 
               cell_numbers: list | tuple, 
               steps: dict[str: list[int]]
               ) -> ArbinBatch:
    """
    Load the entire B6T10/B6T15 dataset.

    Parameters
    ----------
    ``cell_builder`` \: ``CellBuilder``
        ``CellBuilder`` used to access data reading methods.
    ``prepath`` \: ``str``
        path to the folder containing all battery data.
    ``channel_numbers`` \: ``list | tuple``
        Channel numbers to be loaded.
    ``cell_numbers`` \: ``list | tuple``
        Cell number of channels to be loaded. Required to be
        in order correspondin to the channel numbers.
    ``steps`` \: ``dict[str: list[int]]``
        Steps to be loaded.
    """
    # list for holding processed cells
    arbin_cells = []

    # loop over channel numbers
    for channel_idx, channel in enumerate(channel_numbers):
        print(f'Processing channel {channel}')

        # append new cell to cells processed cells list
        arbin_cells.append(ArbinCell(cell_numbers[channel_idx], channel))

        # make subfolder name in raws folder
        # two folders are needed because the tests were modified part way through
        # hence two seperate folders contain the "different" data.
        T10_folder_name = f'B6T10V0_1_2_3_4_9_10_11_12_13_14_15_16/Channel_{channel}/'
        T15_folder_name = f'B6T15V0_1_2_3_4_9_10_11_12_13_14_15_16/Channel_{channel}/'

        # get directory of the current folder
        T10_directory = os.fsencode(prepath+T10_folder_name)
        T15_directory = os.fsencode(prepath+T15_folder_name)

        # sort directory into chronological order for read_B6_csv_data()
        T10_sorted_dir = sorted(os.listdir(T10_directory), key=lambda file: int(os.fsdecode(file).split('.')[1]))
        T15_sorted_dir = sorted(os.listdir(T15_directory), key=lambda file: int(os.fsdecode(file).split('.')[1]))

        # load the files using cellbuilder
        load_files_from_dir(cell_builder, arbin_cells, T10_sorted_dir, prepath, T10_folder_name, steps, channel_idx)
        load_files_from_dir(cell_builder, arbin_cells, T15_sorted_dir, prepath, T15_folder_name, steps, channel_idx)

    # load cells into a batch
    return ArbinBatch(cells=arbin_cells)

def load_B6T11(cell_builder: CellBuilder, 
               prepath: str, 
               channel_numbers: list | tuple, 
               cell_numbers: list | tuple, 
               steps: dict[str: list[int]]
               ) -> ArbinBatch:
    """
    Load the entire B6T11/B6T16 dataset.

    Parameters
    ----------
    ``cell_builder`` \: ``CellBuilder``
        ``CellBuilder`` used to access data reading methods.
    ``prepath`` \: ``str``
        path to the folder containing all battery data.
    ``channel_numbers`` \: ``list | tuple``
        Channel numbers to be loaded.
    ``cell_numbers`` \: ``list | tuple``
        Cell number of channels to be loaded. Required to be
        in order correspondin to the channel numbers.
    ``steps`` \: ``dict[str: list[int]]``
        Steps to be loaded.
    """
    # list for holding processed cells
    arbin_cells = []

    # loop over channel numbers
    for channel_idx, channel in enumerate(channel_numbers):
        print(f'Processing channel {channel}')

        # append new cell to cells processed cells list
        arbin_cells.append(ArbinCell(cell_numbers[channel_idx], channel))

        # make subfolder name in raws folder
        # two folders are needed because the tests were modified part way through
        # hence two seperate folders contain the "different" data.
        T11_folder_name = f'B6T11V0_6_7/Channel_{channel}/'
        T16_folder_name = f'B6T16V0_6_7/Channel_{channel}/'

        # get directory of the current folder
        T11_directory = os.fsencode(prepath+T11_folder_name)
        T16_directory = os.fsencode(prepath+T16_folder_name)

        # sort directory into chronological order for read_B6_csv_data()
        T10_sorted_dir = sorted(os.listdir(T11_directory), key=lambda file: int(os.fsdecode(file).split('.')[1]))
        T15_sorted_dir = sorted(os.listdir(T16_directory), key=lambda file: int(os.fsdecode(file).split('.')[1]))

        # load the files using cellbuilder
        load_files_from_dir(cell_builder, arbin_cells, T10_sorted_dir, prepath, T11_folder_name, steps, channel_idx)
        load_files_from_dir(cell_builder, arbin_cells, T15_sorted_dir, prepath, T16_folder_name, steps, channel_idx)

    # load cells into a batch
    return ArbinBatch(cells=arbin_cells)

def assign_soh(step: int, soh_step: int, nom_cap: float, batch: ArbinBatch) -> None:
    """
    Calculate and assign the SOH to step ``step`` in an attribute named ``soh``
    for all cells in ``batch``.

    Parameters
    ----------
    ``step`` \: ``int``
        Step number to have soh added to
    ``soh_step`` \: ``int``
        Step containing the full discharge data used to calculate the soh.
    ``nom_cap`` \: ``float``
        Nominal capacity of the battery used for the test in Ah.
    ``batch`` \: ``ArbinBatch``
        Batch of cells containing the desired data.
    """
    for cell in batch:
        for cycle in cell:
            try:
                soh = cycle[soh_step][0]['Discharge_Capacity(Ah)'][-1]/nom_cap
            except IndexError:
                # allow for cycles without a soh_step
                soh = -1
            try:
                cycle[step]
            except IndexError:
                # allow for cycles without a step
                continue
            if not cycle[step]:
                continue
            for cycle_step in cycle[step]:
                cycle_step.soh = soh

def assign_soe(step: int, soe_step: int, nom_e_cap: float, batch: ArbinBatch) -> None:
    """
    Calculate and assign the SOE to step ``step`` in an attribute named ``soe``
    for all cells in ``batch``.

    Parameters
    ----------
    ``step`` \: ``int``
        Step number to have soe added to
    ``soe_step`` \: ``int``
        Step containing the full discharge data used to calculate the soe.
    ``nom_e_cap`` \: ``float``
        Nominal energy capacity of the battery used for the test in Wh.
    ``batch`` \: ``ArbinBatch``
        Batch of cells containing the desired data.
    """
    for cell in batch:
        for cycle in cell:
            try:
                soe = cycle[soe_step][0]['Discharge_Energy(Wh)'][-1]/nom_e_cap
            except IndexError:
                # allow for cycles without a soe_step
                soe = -1
            try:
                cycle[step]
            except IndexError:
                # allow for cycles without a step
                continue
            if not cycle[step]:
                continue
            for cycle_step in cycle[step]:
                cycle_step.soe = soe

def filter_by_soh(steps: list[ArbinStep] | dict[int: list[ArbinStep]], 
                  soh_range: float, 
                  soh_lower: int=77, 
                  soh_upper: int=101
                  ) -> dict[int: list[ArbinStep]]:
    """
    Filter steps in ``steps`` such that only the steps with a SOH within the bounds of 
    lower and upper are retained.

    This function requires that all steps in ``steps`` have a soh attribute; see 
    ``assign_soh``.

    Parameters
    ----------
    ``steps`` \: ``list[ArbinStep] | dict[int: list[ArbinStep]]``
        ArbinSteps containing soh attributes to be sorted
    ``soh_range`` \: ``float``
        Width of SOH bins to sort ``steps`` into.
        If ``soh_range`` = 1 then ``steps`` will be sorted into bins
        from 77-78,78-79....,100-101 etc.
    ``soh_lower`` \: ``int``, optional
        Lower SOH bound for sorting.
    ``soh_upper`` \: ``int``
        Upper SOH bound for sorting.
    
    Returns
    -------
    ``dict[int: list[ArbinStep]]``
        Dictionary with key of the lower bound for each range and 
        values of all the steps that match that soh range.
    """
    if isinstance(steps, dict):
        steps = [elem for sublist in steps.values() for elem in sublist]
        
    soh_filtered = {}
    for lower in range(soh_lower, soh_upper, soh_range):
        upper = (lower + soh_range)/100
        lower = lower/100
        soh_filtered[lower] = list(
            filter(lambda x: (lower < x.soh < upper), steps)
            )
    return soh_filtered

def filter_by_temp(step: int, 
                   temp_range: tuple | list, 
                   batch: ArbinBatch
                   ) -> dict[int: list[ArbinStep]]:
    """
    Filter steps in ``batch`` such that only the steps with a temperature within the bounds of 
    ``temp_range`` are retained.

    This function requires that all steps in ``steps`` have a temperature attribute; see 
    ``assign_temp``.

    Parameters
    ----------
    ``step`` \: ``int``
        Step number in batch object that has the temperature attribute.
    ``temp_range`` \: ``tuple | list``
        Element zero is the lower temperature bound and element 1 is the upper temperature bound.
    ``batch`` \: ``ArbinBatch``
        ArbinBatch that contains the test data. See ``primed_data_processing`` package.

    Returns
    -------
    ``dict[int: list[ArbinStep]]``
        Dictionary with key of the cell channel number and values of list of 
        step data within the temperature range.
    """
    temp_filtered = {}
    for cell in batch:
        temp_filtered[cell.channel_number] = []
        for cycle in cell:
            try:
                if round(cycle[step][0].temperature) > temp_range[0] \
                    and round(cycle[step][0].temperature) < temp_range[1]:
                    temp_filtered[cell.channel_number].append(cycle[step][0])
            except:
                # allow for cycles without step step
                pass
    return temp_filtered

def assign_cell_cycle_numbers(batch: ArbinBatch, step_idx: int=None):
    """
    Assign cell and cycle numbers to the steps for ease of access when working with data.

    Parameters
    ----------
    ``batch`` \: ``ArbinBatch``
        ``ArbinBatch`` object containing the test data.
    ``step_idx`` \: ``int``
        Step number to assign cell and cycle numbers too.
    """
    for cell in batch:
        for cycle in cell:
            if step_idx:
                try:
                    cycle[step_idx]
                except IndexError:
                    # allow for cycles with no step
                    continue
                for step in cycle[step_idx]:
                    step.cell_number = cell.cell_number
                    step.channel_number = cell.channel_number
            else:
                for step in cycle:
                    step.cell_number = cell.cell_number
                    step.channel_number = cell.channel_number

######################################################## GAMRY ######################################################

def load_B6T10_eis(T10_eis_folder: str, 
                   T15_eis_folder: str, 
                   channel_numbers: list | tuple, 
                   cell_numbers: list | tuple
                   ) -> list[EisCell]:
    """
    Load eis data from batch B6 and test 10 and 15 into a data object.

    Parameters
    ----------
    ``T10_eis_folder`` \: ``str``
        Absolute path to the folder containing test 10 eis data.
    ``T15_eis_folder`` \: ``str``
        Absolute path to the folder containing test 15 eis data.
    ``channel_numbers`` \: ``list | tuple``
        Iterable containing the channel numbers to be read into the data object.
    ``cell_numbers`` \: ``list | tuple``
        Iterable containing the cell numbers to be read into the data object.

    Returns
    -------
    ``list[EisCell]``
        List of ``EisCell`` objects for each cell in the folder directories.
    """
    eis_cells = []
    for channel_idx, channel in enumerate(channel_numbers):
        cycle = 1
        eis_cycles = []
        while cycle <= 46:
            eis_sweep = EisSweep(f'eis cycle{cycle}', 0.5, 14)

            try:
                if cycle <= 23:
                    file_prepath = T10_eis_folder
                if cycle > 23:
                    file_prepath = T15_eis_folder
                if cycle < 10 and channel < 10:
                    eis_sweep.read_DTA_file(file_prepath + f'B6T10V0_Chan00{channel}_Cycle00{cycle}_Step014.DTA')
                elif cycle < 10 and channel < 100:
                    eis_sweep.read_DTA_file(file_prepath + f'B6T10V0_Chan0{channel}_Cycle00{cycle}_Step014.DTA')
                elif cycle < 100 and channel < 10:
                    eis_sweep.read_DTA_file(file_prepath + f'B6T10V0_Chan00{channel}_Cycle0{cycle}_Step014.DTA')
                elif cycle < 100 and channel < 100:
                    eis_sweep.read_DTA_file(file_prepath + f'B6T10V0_Chan0{channel}_Cycle0{cycle}_Step014.DTA')
                else:
                    print('Cycle number greater than 100!')

                eis_cycles.append(EisCycle(cycle, [eis_sweep], f'cycle_object_{cycle}'))
                cycle += 2

            except FileNotFoundError:
                cycle += 2
                print(f'Warning... File B6T10V0_Chan00{channel}_Cycle00{cycle}_Step014.DTA doesn\'t exist!')
            
        eis_cells.append(EisCell(
            name=f'eis step for channel{channel}', 
            eis_cycles=eis_cycles, 
            cell_number=cell_numbers[channel_idx], 
            channel_number=channel)
            )
        eis_cycles = []
    
    return eis_cells

def assign_temp(eis_step: int, 
                temp_step: int, 
                batch: ArbinBatch,
                default_step_idx: int=-1,
                default_temp_idx: int=-1
                ) -> None:
    """
    Assign cell temperature to the eis step from the temperature of ``temp_step``.

    Parameters
    ----------
    ``eis_step`` \: ``int``
        Step number of the eis step
    ``temp_step`` \: ``int``
        Step number of the step with temperature data.
    ``batch`` \: ``ArbinBatch``
        ``ArbinBatch`` object with the test data.
    ``default_step_idx`` \: ``int``, optional
        Index of the list of steps within the cycle. Default is -1.
    ``default_temp_idx`` \: ``int``, optional
        Index of the list of temperatures within the step. Default is -1.
    """
    for cell in batch:
        for cycle in cell:
            try:
                # assign the last temperature of the step before EIS to the eis step
                if len(cycle[temp_step]) > 1:
                    print(f"WARNING: More than one temperature step in cycle {cycle.cycle_index}. The last step will be taken by default.")
                temperature = cycle[temp_step][default_step_idx]['Battery_Temperature(C)'][default_temp_idx]
            except IndexError:
                temperature = None
            try:
                cycle[eis_step]
            except IndexError:
                continue
            for step in cycle[eis_step]:
                step.temperature = temperature

def get_first_quadrant_data(step: int, batch: ArbinBatch) -> None:
    """
    Make a dictionary containing only the first quadrant eis data.
    See ``ArbinStep.make_first_quadrant_dict()``.

    Parameters
    ----------
    ``step`` \: ``int``
        Step number to get first quadrant data from.
    ``batch`` \: ``ArbinBatch``
        ``ArbinBatch`` object containing the test data.
    """
    for step in batch[:,:,step]:
        step.make_first_quadrant_dict()
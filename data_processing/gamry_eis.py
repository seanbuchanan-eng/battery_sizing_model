# Sean Buchanan
# This module contains classes that can be used to process and store EIS data
# for ease of use.

import numpy as np
import pandas as pd

class EisSweep:
    """
    Represents a single EIS sweep. Made specifically to represent EIS from 
    the B6 batch of data at PRIMED.

    Attributes
    ----------
    name : str
        Name of the eis sweep.
    soc : float
        State-of-charge of battery when eis measurement was taken. 0 < soc <=1
    pt : list
        Point numbers of the measurements in a Gamry DTA file.
    time : list
        Timesteps of the measurements in a Gamry DTA file in seconds.
    freq : list
        Frequency of the measurements in a Gamry DTA file in Hz.
    z_real : list
        Real portion of the batteries measured impedance frequency response in Ohms.
    z_imag : list
        Imaginary portion of the batteries measured impedance frequency response in Ohms.
    z_sig : list
        Honestly not sure what this is other than a column in a Gamry DTA file in Volts.
    z_mod : list
        Magnitude of the batteries measured impedance frequency response in Ohms.
    z_phase : list
        Phase of the batteries measured impedance frequency response in degrees.
    idc : list
        Measured DC current in the Gamry DTA file in Amps.
    vdc : list
        Measured DC voltage in the Gamry DTA file in Volts.
    ie_reange : list
        Not sure. For now it is just a column in the Gamry DTA file.
        
    Methods
    -------
    """

    def __init__(self, name, soc):
        """
        Parameters
        ----------
        name : str
            Name of the EisSweep
        soc : float
            State-of-charge of battery when eis measurement was taken. 0 < soc <=1
        """
        self.name = name
        self.soc = soc
        self.pt = []
        self.time = []
        self.freq = []
        self.z_real = []
        self.z_imag = []
        self.z_sig = []
        self.z_mod = []
        self.z_phase = []
        self.idc = []
        self.vdc = []
        self.ie_range = []
        self._headers = []
        self._data_already_read = False

    def __str__(self):
        pass

    def read_DTA_file(self, file_path) -> None:
        """
        Read a Gamry DTA output file and store its information.

        Parameters
        ----------
        file_path : str
            Path to the data file.
        """
        if self._data_already_read:
            print('Data has already been read for this object!')
            return

        with open(file_path) as file:
            read = False
            for line in file:
                if line.startswith('	Pt'):
                    self._headers = line.split()
                elif line.startswith('	#'):
                    units = line.split()
                    for idx, header in enumerate(self._headers):
                        self._headers[idx] = header + ' (' + units[idx] + ')'
                elif line.startswith('	0'):
                   read = True
                if read == True:
                    data = line.split()
                    self.pt.append(int(data[0]))
                    self.time.append(float(data[1]))
                    self.freq.append(float(data[2]))
                    self.z_real.append(float(data[3]))
                    self.z_imag.append(float(data[4]))
                    self.z_sig.append(float(data[5]))
                    self.z_mod.append(float(data[6]))
                    self.z_phase.append(float(data[7]))
                    self.idc.append(float(data[8]))
                    self.vdc.append(float(data[9]))
                    self.ie_range.append(int(data[10]))

        self._data_already_read = True

    def get_data_as_array(self) -> np.array:
        """
        Returns
        -------
        numpy array
            Nx11 array where the columns are pt, time, freq, zreal, zimag, zsig, zmod,
            zphase, idc, vdc, ierange
        """
        return np.array(
            [self.pt, self.time, self.freq, self.z_real, self.z_imag,
            self.z_sig, self.z_mod, self.z_phase, self.idc, self.vdc, 
            self.ie_range]).T

    def get_data_as_dataframe(self) -> pd.DataFrame:
        """
        Returns
        -------
        pandas DataFrame
            Nx11 table where column headers are pt, time, freq, zreal, zimag, zsig,
            zmod, zphase, idc, vdc, ierange
        """
        return pd.DataFrame(self.get_data_as_array(), columns=self._headers)
        

class EisCycle:
    """
    Represents all of the eis sweeps in a single test cycle.
    """
    
    def __init__(self, cycle_number, sweeps=[], name=None) -> None:
        """
        Parameters
        ----------
        cycle_number : int
            The cycle number in the test.
        sweeps : list, optional
            EisSweep objects (default is empty). It is recommended to put them in chronological order.
        name : str, optional
            Name of the eis cycle (defaule is None)
        """

        self.cycle_number = cycle_number
        self.sweeps = sweeps
        self.name = name
        self.soh = None
        self.battery_temperature = None

    def __str__(self) -> str:
        pass

    def __iter__(self):
        self.iter_index = 0
        return self
        
    def __next__(self) -> EisSweep:
        if self.iter_index < len(self.sweeps):
            eis_sweep = self.sweeps[self.iter_index]
            self.iter_index += 1
            return eis_sweep
        else:
            raise StopIteration

    def __getitem__(self, idx: int) -> EisSweep:
        return self.sweeps[idx]

    def add_sweep(self, eis_sweep) -> None:
        """
        Add an EisSweep object to the cycle.

        Parameters
        ----------
        eis_sweep : EisSweep
            EisSweep object to be added to the cycle.
        """
        self.sweeps.append(eis_sweep)

    def get_data_as_array(self) -> np.array:
        """
        Returns
        -------
        numpy array
            All EisSweeps in the cycle as a m x n x 11 array where m is the number of sweeps, n is the 
            number of measurements, and 11 is the number of columns.

            Can also be seen as an output of m EisSweep.get_data_as_array().
        """
        output_list = []
        for eis_sweep in self.sweeps:
            output_list.append(eis_sweep.get_data_as_array())
        
        return np.array(output_list)

class EisCell:
    """
    Represents all of the the EIS measurements over a single cells lifetime.
    """

    def __init__(self, name, eis_cycles=[], cell_number=None, channel_number=None) -> None:
        """
        Parameters
        ----------
        name : str
            Name of the cell.
        eis_cycles : list, optional
            EisCycle objects representing all of the sweeps measured in a given cycle 9def.
        cell_number : int, optional
            The cell number from the test.
        channel_number : int, optional
            The channel number from the test.
        """
        self.name = name
        self.eis_cycles = eis_cycles
        self.cell_number = cell_number
        self.channel_number = channel_number

    def __iter__(self):
        self.iter_index = 0
        return self
        
    def __next__(self) -> EisSweep:
        if self.iter_index < len(self.eis_cycles):
            eis_cycle = self.eis_cycles[self.iter_index]
            self.iter_index += 1
            return eis_cycle
        else:
            raise StopIteration

    def __getitem__(self, idx: int) -> EisCycle:
        return self.eis_cycles[idx]

    def add_cycle(self, eis_cycle) -> None:
        """
        Parameters
        ----------
        eis_cycle : EisCycle
            EIS cycle to be added to the cell. It's recommended to add them in chronological order.
        """
        self.eis_cycles.append(eis_cycle)

    def get_data_as_array(self) -> np.array:
        """
        Returns
        -------
        numpy array
            All EisCycles that the cell has experienced. Shape is z x m x n x 11 where z is the number
            of cycles, m is the number of sweeps, and n is the number of data points in the sweep. 11 is
            the number of measured parameters in the test.
        """
        output_list = []
        for eis_cycle in self.eis_cycles:
            output_list.append(eis_cycle.get_data_as_array())
        
        return np.array(output_list)
        

if __name__ == '__main__':

    # test get_data_as_array()
    eis_file_path = 'data_processing/testing_data/B6T10V0_Chan001_Cycle001_Step014.DTA'
    test_sweep = EisSweep('test_sweep')
    test_sweep.read_DTA_file(eis_file_path)

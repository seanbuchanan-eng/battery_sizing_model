"""
This script calculates aggregate parameters used in the state of health 
estimation model. The aggregate parameters are necessary to reduce the amount
of data used, hence, the time needed for fitting of the state of health estimation
model.

It relies on ``save_batch_to_disk.py`` to be run first so that the data objects are
available to be loaded into this script for processing.
"""

import os
os.chdir('..')
print(os.getcwd())
import pickle
import numpy as np
import primed_data_processing.primed_utils as primed_utils
from primed_data_processing.arbin_cycler import ArbinBatch

def open_data_file(cell_number, test_number):
    with open(f"data_objects/B6T{test_number}C{cell_number}_active_steps.pkl", "rb") as f:
        batch = ArbinBatch(cells=pickle.load(f))
    return batch

class EFC:

    def __init__(self) -> None:
        self.soc_chg = 0
        self.soc_dchg = 0
        self.cycle_chg = 0
        self.cycle_dchg = 0
    
    def increment_efc(self, step, delta_t=5, nom_cap=4.0, charge_eff=1, discharge_eff=1):
        for idx, c in enumerate(step["Current(A)"]):
            if idx == 0:
                delta_t = step["Step_Time(s)"][idx]
            else:
                delta_t = step["Step_Time(s)"][idx] - step["Step_Time(s)"][idx-1]
            if c < 0:
                self.soc_dchg += np.abs(c) * delta_t / (3600 * nom_cap) * discharge_eff
                if self.soc_dchg >= 1:
                    self.cycle_dchg += 1
                    self.soc_dchg = 0
            else:
                self.soc_chg += c * delta_t / (3600 * nom_cap) * charge_eff
                if self.soc_chg >= 1:
                    self.cycle_chg += 1
                    self.soc_chg = 0

class DOD:

    def __init__(self, soc_init, ref_soc=0.5) -> None:
        self.dod = []
        self.soc = soc_init
        self.delta_soc = 0
        self.soc_sign = 1
        self.ref_soc = ref_soc

    def increment_dod(self, step, delta_t=5, nom_cap=4.0, charge_eff=1, discharge_eff=1):
            init_soc = self.soc
            self.calc_soc(step, delta_t, nom_cap, charge_eff, discharge_eff)
            self.delta_soc = self.soc - init_soc

            if (np.sign(self.delta_soc) != np.sign(self.soc_sign)  
                and np.sign(self.soc_sign) != 0):

                self.soc_sign *= -1
                self.dod.append(np.abs(1 - (1 - self.soc) / self.ref_soc))

    def calc_soc(self, step, delta_t, nom_cap, charge_eff, discharge_eff):
        for idx, c in enumerate(step["Current(A)"]):
            if idx == 0:
                delta_t = step["Step_Time(s)"][idx]
            else:
                delta_t = step["Step_Time(s)"][idx] - step["Step_Time(s)"][idx-1]
            if c < 0:
                self.soc = self.soc + (c * delta_t / (3600 * nom_cap)) / discharge_eff
            else:
                self.soc = self.soc + (c * delta_t / (3600 * nom_cap)) * charge_eff
            if self.soc > 1:
                self.soc = 1.0
            elif self.soc < 0:
                self.soc = 0
                
    def get_dod(self):
        if self.dod[0] == 0:
            return self.dod[1:]
        return self.dod

def calc_fit_factors(batch):
    efc = {}
    mean_c_rate = {}
    mean_dod = {}
    time = {}
    soc_init = 0.5
    soh = get_soh_dict(batch)

    for cycle in batch[:,:]:
        try:
            ref_soh = soh[cycle.cycle_index]/2
        except:
            pass
        efc_cycle = EFC()
        dod_cycle = DOD(soc_init=soc_init, ref_soc=ref_soh)
        c_rate = np.array([])
        cycle_time = 0

        for idx, step in enumerate(cycle):
            if step.step_index == 4:
                dod_cycle.soc = 0
            delta_t = step["Test_Time(s)"][-1] - step["Test_Time(s)"][0]
            if delta_t <= 5:
                continue
            efc_cycle.increment_efc(step)
            dod_cycle.increment_dod(step)
            if np.mean(np.abs(step["Current(A)"])/4.0) < 0.28:
                pass
            else:
                c_rate = np.append(c_rate, np.abs(step["Current(A)"])/4.0)
            cycle_time = step["Test_Time(s)"][-1]/3600

        if not dod_cycle.get_dod():
            dod = 0
        else:
            dod = dod_cycle.get_dod()

        soc_init = dod_cycle.soc

        mean_c_rate[cycle.cycle_index] = np.mean(c_rate)
        mean_dod[cycle.cycle_index] = np.mean(dod)
        efc[cycle.cycle_index] = (efc_cycle.cycle_chg + efc_cycle.cycle_dchg) / 2
        time[cycle.cycle_index] = cycle_time

    return efc, mean_c_rate, mean_dod, time

def get_soh_dict(batch):
    # assign SOH to step 10 in batch
    # must remember that SOH is calculated after
    # the step has gone through all of it's cycling
    primed_utils.assign_soh(10,10,4,batch)

    soh = {}
    for cycle in batch[:,:]:
        try:
            step_10 = cycle[10]
            soh[cycle.cycle_index] = step_10[0].soh
        except IndexError:
            pass
    return soh

if __name__ == '__main__':
    cell_numbers = (9,10,11,12,1,2,3,4,5,6,7,8)
    test_number = 1015
    # cell_numbers = (10,)
    # cell_numbers = (13,16)
    # test_number = 1217
    # cell_numbers = (14,15)
    # test_number = 1116
    for cell_number in cell_numbers:
        batch = open_data_file(cell_number, test_number)

        efc, mean_c, dod, time = calc_fit_factors(batch)
        soh = get_soh_dict(batch)

        with open(f"aggregate_params/C{cell_number}.npy", "wb") as f:
            np.save(f, np.array([efc, mean_c, dod, time, soh]))

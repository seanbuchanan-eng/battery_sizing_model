"""
This module contains all of the classes that represent distinct blocks
within a battery simulation.

@author Sean Buchanan
"""

import pickle
import numpy as np

class ECM:
    """
    Represents an equivalent circuit model.

    This class is used for simulating the circuit given
    circuit parameters. Other objects handle the determination
    of fitted parameters from data.

    Currently only supports multi-RC Randles circuits.
    """

    def __init__(self, num_rc_pairs: int, R0: float, R: list[float], C: list[float], timestep: int) -> None:

        self.r0 = R0
        self.R = np.array(R)
        self.C = np.array(C)
        self.A_rc = np.zeros((len(R), len(R)))
        self.B_rc = np.zeros((len(R), 1))
        self.I_rk = np.ones((len(R), 1))
        for idx in range(num_rc_pairs):
            self.A_rc[idx,idx] = np.exp(-timestep / (R[idx] * C[idx]))
            self.B_rc[idx,0] = 1 - np.exp(-timestep / (R[idx] * C[idx]))

    def step(self, I_k: float, ocv: float) -> float:
        v_k = ocv + np.sum(self.R * self.I_rk.ravel()) + self.r0 * I_k
        self.I_rk = self.A_rc @ self.I_rk + self.B_rc * I_k
        return v_k

class OCVTable:
    """
    Represents an OCV lookup table.
    """

    def __init__(self) -> None:
        with open('lookup_tables/charge_ocv.pkl', 'rb') as f:
            self.charge_ocv = pickle.load(f)
        with open('lookup_tables/discharge_ocv.pkl', 'rb') as f:
            self.discharge_ocv = pickle.load(f)
        with open('lookup_tables/soc_interp_points.pkl', 'rb') as f:
            self.x = pickle.load(f)

    def get_charge_ocv(self, soh, soc):
        soh = np.floor(soh * 100) / 100
        ocv = np.interp(soc, self.x, self.charge_ocv[soh])
        return ocv
    
    def get_discharge_ocv(self, soh, soc):
        soh = np.floor(soh * 100) / 100
        ocv = np.interp(soc, self.x, self.discharge_ocv[soh])
        return ocv

class SOCIntegrator:
    """
    Represents and integrator for counting couloumbs.
    """

    def __init__(self, charge_eff: float=1, discharge_eff: float=1) -> None:
        self.charge_eff = charge_eff
        self.discharge_eff = discharge_eff

    def step(self, soc_k, timestep, I_k, nom_cap) -> float:
        """
        Assumes that current is constant over the sampling interval
        """
        if I_k < 0: eff = self.discharge_eff
        else: eff = self.charge_eff
        # if abs(I_k) < 0.01: I_k = 0
        soc_k = soc_k + (I_k * eff * timestep/3600) / nom_cap
        if soc_k >= 1:
            return 1
        elif soc_k <= 0:
            return 0
        return soc_k
    
class DegradationModel:
    """
    Represents a model that predicts the battery SOH.
    """

    def __init__(self, ref_soc=0.5) -> None:
        self.model_params = np.array([4.91441030e-04, -2.98808476e-03,  
                                      3.09008644e-01,  2.80785041e+00,
                                      6.12442315e-01,  1.00000000e-01,  
                                      4.61330493e-01])
        self.dod = []
        self.delta_soc = 0
        self.soc_sign = 1
        self.ref_soc = ref_soc
        self.soc_chg = 0
        self.cycle_chg = 0

    def model(self, efc, c_rate, dod, t, soh_i):
        a = self.model_params[0]
        b = self.model_params[1]
        c = self.model_params[2]
        d = self.model_params[3]
        e = self.model_params[4]
        alpha = self.model_params[5]
        beta = self.model_params[6]
        mean_c, mean_dod = self._get_window_means(c_rate, dod)
        a_mean_c = np.mean(a * mean_c)
        c_mean_dod = np.mean(c * mean_dod)
        return soh_i*np.exp(b * efc * (c_mean_dod ** d + a_mean_c ** e)) - alpha * t ** beta
    
    def _get_window_means(self, c_rate: list[float], dod: list[float]) -> tuple[list[float], list[float]]:
        """
        Assumes that c_rate and dod are the same length.
        """
        if len(c_rate) == 1:
            c_window = np.array([c_rate, c_rate, c_rate, c_rate])
            dod_window = np.array([dod, dod, dod, dod])
        elif len(c_rate) == 2:
            c_window = np.array([c_rate[-2], c_rate[-1], c_rate[-1], c_rate[-1]])
            dod_window = np.array([dod[-2], dod[-1], dod[-1], dod[-1]])
        elif len(c_rate) == 3:
            c_window = np.array([c_rate[-3], c_rate[-2], c_rate[-1], c_rate[-1]])
            dod_window = np.array([dod[-3], dod[-2], dod[-1], dod[-1]])
        else:
            c_window = np.array([c_rate[-4], c_rate[-3], c_rate[-2], c_rate[-1]])
            dod_window = np.array([dod[-4], dod[-3], dod[-2], dod[-1]])
        return c_window, dod_window

    def increment_dod(self, prev_soc, curr_soc):
        self.delta_soc = curr_soc - prev_soc

        if (np.sign(self.delta_soc) != np.sign(self.soc_sign)  
            and np.sign(self.soc_sign) != 0):

            self.soc_sign *= -1
            self.dod.append(np.abs(1 - (1 - curr_soc) / self.ref_soc))

    def get_dod(self):
        if self.dod[0] == 0:
            return self.dod[1:]
        return self.dod
    
    def increment_efc(self, prev_soc, curr_soc):
        self.soc_chg += np.abs(prev_soc - curr_soc)
        if self.soc_chg >= 1:
            self.cycle_chg += 1
            self.soc_chg = 0
    
    def get_efc(self):
        # meant to be (charge cycles + discharge cycles)/2 but
        # charge cycles and discharge cycles have been combined
        # in increment_efc hence the / 2 here.
        return self.cycle_chg / 2

class BatterySimulator:
    """
    Represents a battery simulator.
    """

    def __init__(self, 
                 ecm: ECM, 
                 ocv: OCVTable, 
                 integrator: SOCIntegrator,
                 deg_model: DegradationModel, 
                 soc: float=1, 
                 soh: float=1, 
                 nom_cap: float=4.0,
                 ) -> None:
        self.voltage = ocv.get_discharge_ocv(soh, soc)
        self.soc = soc
        self.soh = soh
        self.soh_i = soh
        self.ecm = ecm
        self.ocv = ocv
        self.integrator = integrator
        self.deg_model = deg_model
        self.nom_cap = nom_cap
        self.voltage_result = []
        self.soc_result = []
        self.soh_result = [soh]

    def run(self, current, timestep) -> None:
        # soh model data
        time = 0
        total_time = 0
        mean_c = []
        mean_dod = []
        prev_soc = 0.5
        c_rate = []

        # event loop
        for I in current:
            # ECM
            if I > 0:
                ocv = self.ocv.get_charge_ocv(self.soh, self.soc)
            else:
                ocv = self.ocv.get_discharge_ocv(self.soh, self.soc)
            self.voltage = self.ecm.step(I, ocv)
            prev_soc = self.soc
            self.soc = self.integrator.step(self.soc, timestep, I, self.nom_cap)
            self.voltage_result.append(self.voltage)
            self.soc_result.append(self.soc)

            #SOH Model
            time += timestep/3600 # time in hours
            # increment mean c-rate
            if np.abs(I) > 0.1: c_rate.append(np.abs(I)/4.0)
            # increment mean DOD
            self.deg_model.increment_dod(prev_soc, self.soc, )
            # increment efc (total cumlative)
            self.deg_model.increment_efc(prev_soc, self.soc)

            if time > 340:
                # update mean c-rate for this time window
                mean_c.append(np.mean(c_rate))
                # update mean DOD for this time window
                mean_dod.append(np.mean(self.deg_model.get_dod()))
                # increment total cumlative time
                total_time += time
                # calculate current soh (need to work in backward looking window)
                self.soh = self.deg_model.model(self.deg_model.get_efc(), mean_c, mean_dod, total_time, self.soh_i*100)/100
                # start the next step
                time = 0
                c_rate = []
                self.deg_model.dod = []
                self.soh_result.append(self.soh)

if __name__ == "__main__":
    ecm = ECM(3, 1,  [1,1,1], [1,1,1], 5)
    print(ecm.A_rc)
    print(ecm.B_rc)
    print(ecm.I_rk)
    print(f"R * I_rk: {np.sum(ecm.R * ecm.I_rk)}")
    print(ecm.A_rc @ ecm.I_rk)
    print(ecm.B_rc * 2)

    v_k = ecm.step(2, 3.75)
    print(v_k)
    print(ecm.I_rk)



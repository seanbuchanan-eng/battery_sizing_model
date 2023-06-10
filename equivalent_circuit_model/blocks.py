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
        return soc_k + (I_k * eff * timestep/3600) / nom_cap

class BatterySimulator:
    """
    Represents a battery simulator.
    """

    def __init__(self, 
                 ecm: ECM, 
                 ocv: OCVTable, 
                 integrator: SOCIntegrator, 
                 soc: float=1, 
                 soh: float=1, 
                 nom_cap: float=4.0,
                 ) -> None:
        self.voltage = ocv.get_discharge_ocv(soh, soc)
        self.soc = soc
        self.soh = soh
        self.ecm = ecm
        self.ocv = ocv
        self.integrator = integrator
        self.nom_cap = nom_cap
        self.voltage_result = []
        self.soc_result = []

    def run(self, current, timestep) -> None:
        # event loop
        for I in current:
            if I > 0:
                ocv = self.ocv.get_charge_ocv(self.soh, self.soc)
            else:
                ocv = self.ocv.get_discharge_ocv(self.soh, self.soc)
            self.voltage = self.ecm.step(I, ocv)
            self.soc = self.integrator.step(self.soc, timestep, I, self.nom_cap)
            self.voltage_result.append(self.voltage)
            self.soc_result.append(self.soc)


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



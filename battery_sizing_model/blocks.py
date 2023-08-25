"""
This module contains all of the classes that represent distinct blocks
within a battery simulation.
"""

import pickle
import numpy as np

class CapacityError(ValueError):
    
    def __init__(self, 
                 *args: object, 
                 message="State of charge is greater than the current state of health or less than 0."
                 ) -> None:
        self.message = message
        super().__init__(self.message, *args)

class ECM:
    """
    Represents an equivalent circuit model.

    This class is used for simulating a circuit given
    circuit parameters.

    Currently only supports multi-RC Randles circuits.

    Attributes
    ----------
    ``r0`` \: ``float``
        Series resistance of the Randles circuit.
    ``R`` \: ``numpy.array``
        Array of resistances in the RC pairs.
    ``C`` \: ``numpy.array``
        Array of capacitances in the RC pairs.
    ``A_rc`` \: ``numpy.array``
        State matrix for the circuit.
    ``B_rc`` \: ``numpy.array``
        Control matrix for the circuit.
    ``I_rk`` \: ``numpy.array``
        Previous current estimate through the RC pairs.
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
        """
        Steps the circuit model forward one timestep.

        Parameters
        ----------
        ``I_k`` \: ``float``
            Current input for the this step.
        ``ocv`` \: ``float``
            Open circuit voltage from the ocv table.

        Returns
        -------
        ``float``
            Predicted battery voltage.
        """
        self.I_rk = self.A_rc @ self.I_rk + self.B_rc * I_k
        v_k = ocv + np.sum(self.R * self.I_rk.ravel()) + self.r0 * I_k
        return v_k

class OCVTable:
    """
    Represents an open circuit voltage lookup table.

    Currently only supports the molicel P42A cell (NCA).

    Attributes
    ----------
    ``charge_ocv`` \: ``dict``
        Charge open circuit voltage indexed by state of health.
    ``discharge_ocv`` \: ``dict``
        Discharge open circuit voltage indexed by state of health.
    ``x`` \: ``array``
        Interpolation points for the lookup table.
    """

    def __init__(self) -> None:
        with open('lookup_tables/charge_ocv.pkl', 'rb') as f:
            self.charge_ocv = pickle.load(f)
        with open('lookup_tables/discharge_ocv.pkl', 'rb') as f:
            self.discharge_ocv = pickle.load(f)
        with open('lookup_tables/soc_interp_points.pkl', 'rb') as f:
            self.x = pickle.load(f)

    def get_charge_ocv(self, soh: float, soc: float) -> float:
        """
        Gets the charge open circuit voltage from a lookup table.

        Parameters
        ----------
        ``soh`` \: ``float``
            Current state of health of the battery.
        ``soc`` \: ``float``
            Current state of charge of the battery.

        Returns
        -------
        ``float``
            Open circuit voltage from a lookup table.
        """
        soh = np.floor(soh * 100) / 100
        # 0.96 because thats all the data we have available for the molicel battery
        if soh > 0.96: soh = 0.96
        ocv = np.interp(soc, self.x, self.charge_ocv[soh])
        return ocv
    
    def get_discharge_ocv(self, soh: float, soc: float) -> float:
        """
        Gets the discharge open circuit voltage from a lookup table.

        Parameters
        ----------
        ``soh`` \: ``float``
            Current state of health of the battery.
        ``soc`` \: ``float``
            Current state of charge of the battery.

        Returns
        -------
        ``float``
            Open circuit voltage from a lookup table.
        """
        soh = np.floor(soh * 100) / 100
        # 0.96 because thats all the data we have available for the molicel battery
        if soh > 0.96: soh = 0.96
        ocv = np.interp(soc, self.x, self.discharge_ocv[soh])
        return ocv

class SOCIntegrator:
    """
    Represents an integrator for counting couloumbs and calculating state of charge.

    Attributes
    ----------
    ``charge_eff`` \: ``float``, optional
        Charge efficiency of the battery [0,1]. Default = 1.
    ``discharge_eff`` \: ``float``, optional
        Discharge efficiency of the battery [0,1]. Default = 1.
    ``validate`` \: ``bool``, optional
        Flag for indicating whether the sim should compensate for the
        varying charge/discharge efficiency in real data when validating.
        Default = False.
    """

    def __init__(self, charge_eff: float=1, discharge_eff: float=1, validate: bool=False) -> None:
        self.charge_eff = charge_eff
        self.discharge_eff = discharge_eff
        self.validate = validate

    def step(self, soc_k: float, timestep: float, I_k: float, nom_cap: float) -> float:
        """
        Integrate the current for this timestep.

        Assumes that current is constant over the sampling interval.

        Parameters
        ----------
        ``soc_k`` \: ``float``
            State of charge of the previous step.
        ``timestep`` \: ``float``
            Duration in time for this step.
        ``I_k`` \: ``float``
            Current for this step.
        ``nom_cap`` \: ``float``
            Nominal capacity of the battery.

        Returns
        -------
        ``float``
            State of charge from this step.
        """
        if I_k < 0: eff = 1/self.discharge_eff
        else: eff = self.charge_eff
        soc_k = soc_k + (I_k * eff * timestep/3600) / nom_cap

        # Used to compensate for the varying charge/discharge efficiency
        # in real data when validating the simulator.
        if self.validate and soc_k >= 1:
            return 1
        elif self.validate and soc_k <= 0:
            return 0
        return soc_k
    
class DegradationModel:
    """
    Represents a model that predicts the battery SOH.

    Attributes
    ----------
    ``model_params`` \: ``numpy.array``
        Fitted coefficients for the model using molicel P42A cell data.
    ``dod`` \: ``list``
        Depth of discharge of every step in the simulation.
    ``delta_soc`` \: ``float``
        Change in state of charge between steps.
    ``soc_sign`` \: ``int``
        1 for positive change in state of charge. -1 for negative change in state of charge.
    ``ref_soc`` \: ``float``
        Reference state of charge for the depth of discharge calculation.
    ``soc_chg`` \: ``float``
        Amount of charge added to the battery for each step.
    ``cycle_chg`` \: ``int``
        Number of equivalent full cycles * 2.
    """

    def __init__(self, ref_soc: float=0.5) -> None:
        # model parameters determined from fitting to data.
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

    def model(self, efc: float, c_rate: list[float], dod: list[float], t: float, soh_i: float) -> float:
        """
        State of health estimation model that estimates the state of health of the battery
        for a given time-period.

        Parameters
        ----------
        ``efc`` \: ``float``
            Equivalent full cycles, defined as the amount of charge needed to fully charge
            or discharge the battery. See ``increment_efc()``.
        ``c_rate`` \: ``list[float]``
            C-rate for every soh estimation step.
        ``dod`` \: ``list[float]``
            Depth of discharge for every soh estimation step. See ``increment_dod()``.
        ``t`` \: ``float``
            Total time in hours that the battery has existed.
        ``soh_i`` \: ``float``
            Beginning of life state of health of the battery.

        Returns
        -------
        ``float``
            Estimeated state of health for the given input parameters.
        """
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
        Calculates the following window for the soh model.
        
        Assumes that c_rate and dod are the same length.

        Parameters
        ----------
         ``c_rate`` \: ``list[float]``
            C-rate for every soh estimation step.
        ``dod`` \: ``list[float]``
            Depth of discharge for every soh estimation step. See ``increment_dod()``.

        Returns 
        -------
        ``tuple[list[float], list[float]]``
            A tuple of the last four c-rates and dods.
        """
        if len(c_rate) == 1:
            c_window = np.array([c_rate[-1], c_rate[-1], c_rate[-1], c_rate[-1]])
            dod_window = np.array([dod[-1], dod[-1], dod[-1], dod[-1]])
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

    def increment_dod(self, prev_soc: float, curr_soc: float) -> None:
        """
        Increments the depth of discharge based on the previoius and current soc.

        DOD is only calculated if the sign of the current changes.

        Parameters
        ----------
        ``prev_soc`` \: ``float``
            Previous steps state of charge.
        ``curr_soc`` \: ``float``
            Current steps state of charge.
        """
        self.delta_soc = curr_soc - prev_soc

        if (np.sign(self.delta_soc) != np.sign(self.soc_sign)  
            and np.sign(self.soc_sign) != 0):

            self.soc_sign *= -1
            self.dod.append(np.abs(1 - (1 - curr_soc) / self.ref_soc))

    def get_dod(self):
        """
        Getter for depth of discharge.
        """
        if self.dod[0] == 0:
            return self.dod[1:]
        return self.dod
    
    def increment_efc(self, prev_soc: float, curr_soc: float) -> None:
        """
        Increments the equivalent full cycle value given the previous and current
        state of charge.

        An equivalent full cycle is counted when the battery has experienced enough 
        cumulative charge to fully charge or discharge the battery.

        Parameters
        ----------
        ``prev_soc`` \: ``float``
            Previous steps state of charge.
        ``curr_soc`` \: ``float``
            Current steps state of charge.
        """
        self.soc_chg += np.abs(prev_soc - curr_soc)
        if self.soc_chg >= 1:
            self.cycle_chg += 1
            self.soc_chg = 0
    
    def get_efc(self):
        """
        Getter for equivalent full cycles.
        """
        # meant to be (charge cycles + discharge cycles)/2 but
        # charge cycles and discharge cycles have been combined
        # in increment_efc hence the / 2 here.
        return self.cycle_chg / 2
    
class Battery:
    """
    Represents a full battery pack of a single chemistry.

    Attributes
    ----------
    ``cell_chemistry`` \: ``str``
        Chemistry of the battery. Currently only supports NCA.
    ``pack_capacity`` \: ``float``
        Capacity of the pack in Watt-hours.
    ``pack_voltage`` \: ``float``, optional
        Rated voltage of the pack, doesn't have any bearing on the
        performance of the simulation. Default = 0.
    ``cell_dishcarge_energy_eff`` \: ``float``, optional
        Energy efficiency of the battery during discharge. Default = 1.
    ``cell_charge_energy_eff`` \: ``float``, optional
        Energy efficiency of the battery during charge. Default = 1.
    ``cell_capacity_ah`` \: ``float``
        Cell capacity in Amp-hours.
    ``cell_capacity_wh`` \: ``float``
        Cell capacity in Watt-hours.
    ``cell_nom_voltage`` \: ``float``
        Nominal voltage of the cell.
    ``cell_min_voltage`` \: ``float``
        Minimum voltage of the cell.
    ``price_per_kwh`` \: ``float``
        Estimated price_per_kwh of the cell chemistry.
    ``capital_cost`` \: ``float``
        Capital cost of the battery pack cells.
    ``num_series`` \: ``int``
        Calculated number of cells in series in the pack. Only valid
        if ``pack_voltage`` > 0.
    ``num_parallel`` \: ``int``
        Calculated number of cells in parallel in the pack.
    """

    def __init__(self, 
                 cell_chemistry: str,
                 pack_capacity: float, #Wh
                 pack_voltage: float=0,
                 cell_discharge_energy_eff: float=1,
                 cell_charge_energy_eff: float=1,
                ) -> None:
        
        self.pack_capacity = pack_capacity
        self.pack_voltage = pack_voltage
        self.cell_discharge_energy_eff = cell_discharge_energy_eff
        self.cell_charge_energy_eff = cell_charge_energy_eff
        if cell_chemistry.upper() == "NCA":
            self.cell_capacity_ah = 4.2 #Ah
            self.cell_capacity_wh = 14.7 #Wh
            self.cell_nom_voltage = 3.6 #V
            self.cell_max_voltage = 4.2 #V
            self.cell_min_voltage = 2.5 #V
            self.price_per_kwh = 150 #USD
            self.capital_cost = self.pack_capacity/1e3 * self.price_per_kwh
            self.cell_chemistry = cell_chemistry

        self.num_series = np.ceil(self.pack_voltage / self.cell_min_voltage)
        self.num_parallel = np.ceil(self.pack_capacity / self.cell_capacity_wh)
        
    def get_current(self, power: float) -> float:
        """
        Calculates current in Amps from nominal voltage and input/output power.
        
        Parameters
        ----------
        ``power`` \: ``float``
            Input power signal for the battery at each step in the simulation.

        Returns
        -------
        ``float``
            Calculated cell current using the cell nominal voltage.
        """
        if power < 0: return power / self.num_parallel / self.cell_nom_voltage / self.cell_discharge_energy_eff

        return power / self.num_parallel / self.cell_nom_voltage * self.cell_charge_energy_eff
 
class BatterySimulator:
    """
    Represents a simulation for a battery.

    Attributes
    ----------
    ``ecm`` \: ``ECM``
        Equivalent circuit model object.
    ``ocv`` \: ``OCVTable``
        Open circuit voltgate lookup table object.
    ``integrator`` \: ``SOCIntegrator``
        State of charge integrator.
    ``deg_model`` \: ``DegradationModel``
        State of health estimation model.
    ``battery`` \: ``Battery``
        Battery object.
    ``soc`` \: ``float``
        State of charge at this step.
    ``soh`` \: ``float``
        State of health at this step.
    ``voltage`` \: ``float``
        Voltage at this step.
    ``soh_i`` \: ``float``
        Beginning of life state of health.
    ``nom_cap`` \: ``float``
        Cell nominal capacity.
    ``voltage_result`` \: ``list[float]``
        Estimated voltage for every step in the simulation.
    ``current_result`` \: ``list[float]``
        Calculated current for every step in the simulation.
    ``soc_result`` \: ``list[float]``
        State of charge calculated at every step in the simulation.
    ``soh_result`` \: ``list[float]``
        State of health calculated every 340 hours of simulation.
    """

    def __init__(self, 
                 ecm: ECM, 
                 ocv: OCVTable, 
                 integrator: SOCIntegrator,
                 deg_model: DegradationModel,
                 battery: Battery,
                 soc: float=1, 
                 soh: float=1, 
                 ) -> None:
        self.voltage = ocv.get_discharge_ocv(soh, soc)
        self.soc = soc
        self.soh = soh
        self.soh_i = soh
        self.ecm = ecm
        self.ocv = ocv
        self.integrator = integrator
        self.deg_model = deg_model
        self.battery = battery
        self.deg_model.ref_soc = self.soh/2
        self.nom_cap = battery.cell_capacity_ah
        self.voltage_result = []
        self.current_result = []
        self.soc_result = []
        self.soh_result = [soh]

    def run(self, power: list[float], timestep: float, validate: bool=False) -> None:
        """
        Simulation loop that estimates battery state of health, state of charge,
        and voltage.

        Parameters
        ----------
        ``power`` \: ``array[float]``
            Entire power signal for the battery with magnitudes on the pack level scale.
            Points are spaced out in time with a distance of ``timestep`` between them.
        ``timestep`` \: ``float``
            Timestep of the simulation.
        ``validate`` \: ``bool``
            Flag for turning off error checking to account for varying charge/discharge
            efficiencies in real life data.
        """
        # soh model data
        time = 0
        total_time = 0
        mean_c = []
        mean_dod = []
        prev_soc = 0.5
        c_rate = []

        # event loop
        for p in power:
            I = self.battery.get_current(p)
            self.current_result.append(I)

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

            if (self.soc > self.soh or self.soc < 0) and not validate:
                raise CapacityError
            
            if self.voltage > self.battery.cell_max_voltage and not validate:
                print(f"WARNING: Maximum voltage exceeded! Voltage = {self.voltage}")
            elif self.voltage < self.battery.cell_min_voltage and not validate:
                print(f"WARNING: Minimum voltage exceeded! Voltage = {self.voltage}")

            #SOH Model
            time += timestep/3600 # time in hours
            # increment mean c-rate
            c_rate.append(np.abs(I)/self.battery.cell_capacity_ah)
            # increment mean DOD
            self.deg_model.increment_dod(prev_soc, self.soc)
            # increment efc (total cumlative)
            self.deg_model.increment_efc(prev_soc, self.soc)

            if time > 340:
                # update mean c-rate for this time window
                mean_c.append(np.mean(c_rate))
                # update mean DOD for this time window
                mean_dod.append(np.mean(self.deg_model.get_dod()))
                # increment total cumlative time
                total_time += time
                # calculate current soh
                self.soh = self.deg_model.model(self.deg_model.get_efc(), mean_c, mean_dod, total_time, self.soh_i*100)/100
                self.deg_model.ref_soc = self.soh/2
                # start the next step
                time = 0
                c_rate = []
                self.deg_model.dod = []
                self.soh_result.append(self.soh)
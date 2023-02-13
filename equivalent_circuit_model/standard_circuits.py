import numpy as np
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt

class RandlesWarbug:
    """
    This class represents the impedance of a randles circuit with a Warbug
    element in series with the parallel resistor.
            -----C1-----
            |           |
    ---R1---            -----------
            |           |
            ---R2--W0---

    The equation for the circuit impedance is given by:
    R1 + 1 / (1j*omega*C + 1 / (R2 + Y0/np.sqrt(omega)*(1 - 1j)))

    Attributes
	----------
    R1 : float
        Resistance of the R1 resistor in Ohms.
    R2 : float
        Resistance of the R2 resistor in Ohms.
    C : float
        Capacitance of the capacitor in Farads.
    Y0 : float
        Warbug coefficient in Ohm*sec**-0.5.

	Methods
	-------
    """

    def __init__(self, parameters=[1, 1, 1, 1]) -> None:
        """
        Parameters
        ----------
        parameters : list
            Values for the elements in the circuit in the order [R1, R2, C1, W0].
            Resistance in Ohms, capacitance in Farads, warbug in Ohm*sec**-0.5.
            Default [1,1,1,1].
        """
        self.R1 = parameters[0]
        self.R2 = parameters[1]
        self.C1 = parameters[2]
        self.W0 = parameters[3]

    def __str__(self) -> str:
        pass

    def impedance_model(self, frequency, R1, R2, C, Y0) -> complex:
        """
        Calculate the frequency response of the circuit for given frequencies.

        Parameters
        ----------
        frequency : numpy array
            Frequencies to evaluate the impedance at in Hz.
        R1 : float
            Resistance of the R1 resistor in Ohms.
        R2 : float
            Resistance of the R2 resistor in Ohms.
        C : float
            Capacitance of the capacitor in Farads.
        Y0 : float
            Warbug coefficient in Ohm*sec**-0.5.

        Returns
        -------
        complex array
            Circuit response to every frequency in frequency as an array of complex numbers.
        """
        omega = frequency*2*np.pi
        return R1 + 1 / (1j*omega*C + 1 / (R2 + Y0/np.sqrt(omega)*(1 - 1j)))

    def fit(self, frequency, z_real, z_imag) -> tuple:
        """
        Fit the parameters in the circuit from measured electrochemical impedance
        spectroscopy data using non-linear least squares (marquardt-levenberg).
        
        Using this method will update the original parameters of the model to the fitted
        parameters.

        Parameters
        ----------
        frequency : numpy array
            Frequency that each z_real and z_imag point is measured at in Hz.
        z_real : numpy array
            Real part of the measured impedance data in Ohms.
        z_imag : numpy array
            Imaginary part of the measured impedance data in Ohms.

        Returns
        -------
        tuple
            First element is the fitted parameters as a list. Second element is the
            covariance of those parameters.
        """

        z_both = np.hstack([z_real, z_imag])
        poptboth, pcovboth = curve_fit(
                                        self._concat_real_imag_func,
                                        np.hstack([frequency, frequency]),
                                        z_both,
                                        method='lm'
                                        )

        self.R1 = poptboth[0]
        self.R2 = poptboth[1]
        self.C1 = poptboth[2]
        self.W0 = poptboth[3]

        return (poptboth, pcovboth)

    def _concat_real_imag_func(self, frequency, R1, R2, C, Y0) -> np.array:
        """
        Concatenate the real and imaginary portions of the `impedance model` response together.

        This allows for the `scipy.optimize.curve_fit` function to fit what is effectively three
        dimensional data.

        Parameters
        ----------
        frequency : numpy array
            Frequencies to evaluate the impedance at in Hz.
        R1 : float
            Resistance of the R1 resistor in Ohms.
        R2 : float
            Resistance of the R2 resistor in Ohms.
        C : float
            Capacitance of the capacitor in Farads.
        Y0 : float
            Warbug coefficient in Ohm*sec**-0.5.

        Returns
        -------
        numpy array
            Concatenated real and imaginary portions in dimensions in Ohms.
        """
        N = len(frequency)
        x_real = frequency[:N//2]
        x_imag = frequency[N//2:]
        y_real = np.real(self.impedance_model(x_real, R1, R2, C, Y0))
        y_imag = np.imag(self.impedance_model(x_imag, R1, R2, C, Y0))
        return np.hstack([y_real, y_imag])

    # Plotting functions

    def plot_nyquist(self, frequency, z_real_measured=None, z_imag_measured=None) -> None:
        """
        Display the nyquist plot of the circuit for the given frequencies.

        Parameters
        ----------
        frequency : numpy array
            The frequencies in Hz the model is evaluated and plotted at.
        z_real_measured : numpy array
            The measured real portion of the impedance in Ohms
        z_imag_measured : numpy array
            The measured imaginary portion of th impedance in Ohms
        """
        z_real = self.impedance_model(frequency, self.R1, self.R2, self.C1, self.W0).real
        z_imag = self.impedance_model(frequency, self.R1, self.R2, self.C1, self.W0).imag

        plt.plot(z_real, -1*z_imag, label='estimated', marker='o')

        if (z_real_measured is not None) and (z_imag_measured is not None):
            plt.plot(z_real_measured, -1*z_imag_measured, label='measured', marker='*')

        plt.xlabel('Z_real [ohm]')
        plt.ylabel('Z_imag [ohm]')
        plt.legend()
        plt.show();

    def plot_real(self, frequency, z_real_measured=None) -> None:
        """
        Display the real impedance plot of the circuit for the given frequencies.

        Parameters
        ----------
        frequency : numpy array
            The frequencies in Hz the model is evaluated and plotted at.
        z_real_measured : numpy array
            The measured real portion of the impedance in Ohms
        """
        z_real = self.impedance_model(frequency, self.R1, self.R2, self.C1, self.W0).real

        plt.plot(frequency, z_real, label='estimated', marker='o')

        if (z_real_measured is not None):
            plt.plot(frequency, z_real_measured, label='measured', marker='*')

        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Z_real [ohm]')
        plt.legend()
        plt.show();

    def plot_imag(self, frequency, z_imag_measured=None) -> None:
        """
        Display the real impedance plot of the circuit for the given frequencies.

        Parameters
        ----------
        frequency : numpy array
            The frequencies in Hz the model is evaluated and plotted at.
        z_imag_measured : numpy array
            The measured imaginary portion of th impedance in Ohms
        """
        z_imag = self.impedance_model(frequency, self.R1, self.R2, self.C1, self.W0).imag

        plt.plot(frequency, z_imag, label='estimated', marker='o')

        if (z_imag_measured is not None):
            plt.plot(frequency, z_imag_measured, label='measured', marker='*')

        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Z_imaginary [ohm]')
        plt.legend()
        plt.show();
        
import math 
import numpy as np

# Planck constant (Js).
h = 6.62607015e-34
# Boltzmann constant (J/K).
k = 1.380649e-23
# Speed of light (m/s).
c = 299792458

def blackbody_spectrum_wavel(T, wavel_start = 360, wavel_end = 830, step = 1):
    """Compute the black-body spectrum.

    @param T:           The temperature (K)
    @param wavel_start: The start wavelength (nm)
    @param wavel_end:   The end wavelength (nm)
    @return wavelength, spectrum: The wavelengths and spectrum (W/m^2/nm)."""

    # Conversion factor to convert wavelengths from nm to m.
    factor_nm_to_m = 1.0e-9

    wavelength_nm = np.arange(wavel_start, wavel_end + step, step)
    wavelength_m = factor_nm_to_m * wavelength_nm

    # Note that the output units are W/m^2/nm.
    left = 2 * h * (c * c) / (wavelength_m ** 4) / wavelength_nm
    right = 1 / (np.exp(h * c / (wavelength_m * k * T)) - 1)

    spectrum = left * right 

    return wavelength_nm, spectrum
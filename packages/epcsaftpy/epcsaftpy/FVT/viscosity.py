import numpy as np

# Constants
kb = 1.3806488e-23  # [J/K] Boltzman's constant
Na = 6.022142e23  # [mol-1] Avogadro's Number
R = Na * kb    # [J mol-1 K-1] Ideal gas constant


def viscosity_pure(rho, T, eos, Xass0 = None):
    """
    viscosity for pure components using free volume theory
    (https://www.sciencedirect.com/science/article/pii/S0378381219303024)

    Parameters
    ----------
    rho : float
        molar density [mol/m3]
    T : float
        temperature [K]
    eos : object
        an EoS
    Xass0: array, optional
        Initial guess for the calculation of fraction of non-bonded sites

    Returns
    -------
    ten : float
        viscosity [mPa * s]
    """
    Lv, alpha, B = eos.viscosity_parameters
    Lv *= 1e-10

    if isinstance(T, (int, np.int32, np.int64)):
        T = T * 1.0
        
    P = eos.pressure(rho, T)
    RT = R *  T

    visc = (P + 1e-3 * alpha * rho**2 * eos.Mw) * Lv 
    visc *= np.sqrt(1e-3 * eos.Mw/(3*R*T))
    visc *= 1e3 * np.exp(B*((P + 1e-3 * alpha * rho**2 * eos.Mw) / (rho * RT))**1.5)
    return visc 


def viscosity_mix(x, rho, T, eos, Xass0 = None):
    """
    viscosity for mixtures using free volume theory
    (https://www.sciencedirect.com/science/article/pii/S0378381219303024)

    Parameters
    ----------
    x: array_like
        molar fraction array
    rho : float
        molar density [mol/m3]
    T : float
        temperature [K]
    eos : object
        an EoS
    Xass0: array, optional
        Initial guess for the calculation of fraction of non-bonded sites

    Returns
    -------
    ten : float
        viscosity [mPa * s]
    """
    Lvi, alphai, Bi = eos.viscosity_parameters.T
    
    Lv = np.dot(x, Lvi)*1e-10
    alpha = np.dot(x, alphai)
    B = np.dot(x, Bi)

    P = eos.pressure(x, rho, T)
    RT = R *  T
    Mw = np.dot(eos.Mw, x)

    visc = (P + 1e-3 * alpha * rho**2 * Mw) * Lv 
    visc *= np.sqrt(1e-3 * Mw/(3*R*T))
    visc *= 1e3 * np.exp(B*((P + 1e-3 * alpha * rho**2 * Mw) / (rho * RT))**1.5)
    return visc
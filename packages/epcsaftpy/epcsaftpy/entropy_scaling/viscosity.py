import numpy as np

# Constants
kb = 1.3806488e-23  # [J/K] Boltzman's constant
Na = 6.022142e23  # [mol-1] Avogadro's Number
R = Na * kb    # [J mol-1 K-1] Ideal gas constant


A22, B22, C22, D22, E22, = 1.16145, 0.14874, 0.52487, 0.7732, 2.16178 
F22, R22, S22, W22, P22  = 2.43787, -6.435e-4, 18.0323, -0.7683, 7.27371

def viscosity_pure(rho, T, eos, Xass0 = None):
    """
    viscosity for pure components using an excess entropy scaling approach

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
    if isinstance(T, (int, np.int32, np.int64)):
        T = T * 1.0
        
    Tad = T/(eos.eps / kb)
    Omega22 = A22/Tad**B22 + C22/np.exp(D22 * Tad) + E22 / np.exp(F22 * Tad)
    Omega22 += R22 * Tad**B22 * np.sin(S22 * Tad**W22 - P22)

    visc_ref = np.sqrt(eos.Mw * T * kb / (1000 * Na * np.pi))
    visc_ref *= 5/16
    visc_ref /= eos.sigma**2 * Omega22

    s_star = eos.EntropyR_NrhoT(rho, T, Xass0) / (R * eos.ms)
    svec = np.array([1, s_star, s_star**2, s_star**3])
    lnvisc  = np.dot(svec, eos.viscosity_parameters)
    visc = np.exp(lnvisc) * visc_ref * 1e3
    return visc 


def viscosity_mix(x, rho, T, eos, Xass0 = None):
    """
    viscosity for mixtures using an excess entropy scaling approach

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
    Tad = T/(eos.eps / kb)
    Omega22 = A22/Tad**B22 + C22/np.exp(D22 * Tad) + E22 / np.exp(F22 * Tad)
    Omega22 += R22 * Tad**B22 * np.sin(S22 * Tad**W22 - P22)

    visc_ref = np.sqrt(eos.Mw * T * kb / (1000 * Na * np.pi))
    visc_ref *= 5/16
    visc_ref /= eos.sigma**2 * Omega22

    aux1 = np.outer(visc_ref, 1/visc_ref)
    aux2 = np.outer(1/eos.Mw, eos.Mw)
    phij = aux1**0.5 * aux2**0.25 
    phij += 1
    phij **=2
    phij /= (8*(1 + 1/aux2))**0.5  

    visc_mixref = np.sum(x * visc_ref/np.sum(x * phij, axis = 1))

    Ai, Bi, Ci, Di = eos.viscosity_parameters.T
    mxi = x * eos.ms
    mx = np.sum(mxi)
    auxi = mxi/mx
    s_star = eos.EntropyR_NrhoT(x, rho, T, Xass0) / (R * mx)

    lnvisc_mix = np.dot(x, Ai)
    lnvisc_mix += np.dot(auxi, Bi) * s_star
    lnvisc_mix += np.dot(auxi, Ci) * s_star**2
    lnvisc_mix += np.dot(auxi, Di) * s_star**3

    visc = np.exp(lnvisc_mix) * visc_mixref * 1e3
    return visc
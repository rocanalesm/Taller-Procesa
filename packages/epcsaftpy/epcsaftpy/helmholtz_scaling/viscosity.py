import numpy as np

# Constants
kb = 1.3806488e-23  # [J/K] Boltzman's constant
Na = 6.022142e23  # [mol-1] Avogadro's Number
R = Na * kb    # [J mol-1 K-1] Ideal gas constant


A22, B22, C22, D22, E22, = 1.16145, 0.14874, 0.52487, 0.7732, 2.16178 
F22, R22, S22, W22, P22  = 2.43787, -6.435e-4, 18.0323, -0.7683, 7.27371

def viscosity_pure(rho, T, eos, Xass0 = None):
    """
    viscosity for pure components using an helmholtz scaling approach
    Ind. Eng. Chem. Res. 2021, 60, 9231-9245
    https://doi.org/10.1021/acs.iecr.1c00837


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
    tetha, aa, bb, cc, dd = eos.viscosity_parameters
    if isinstance(T, (int, np.int32, np.int64)):
        T = T * 1.0
        
    Tad = T/(eos.eps / kb)
    Omega22 = A22/Tad**B22 + C22/np.exp(D22 * Tad) + E22 / np.exp(F22 * Tad)
    Omega22 += R22 * Tad**B22 * np.sin(S22 * Tad**W22 - P22)

    visc_ref = np.sqrt(eos.Mw * T * kb / (1000 * Na * np.pi ))
    visc_ref *= 5/16
    visc_ref /= eos.sigma**2 * Omega22

    rhom = Na * rho
    beta, eps_beta, d, dia3, mes3, m2e2s3, deta_drho, Fab, betaxpmu2 = eos.temperature_aux(T)
    eta = deta_drho * rhom
    aHS = ((4 - 3*eta)*eta)/(eta - 1)**2
    

    ares, _ = eos.ares(rhom, T, Xass0 = Xass0)
    aaux = (ares - eos.ms * aHS)
    aaux *= Tad**tetha  

    lnvisc  = aa + bb*aaux
    lnvisc += cc * aaux**2 + dd * aaux**3
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


    # Parameters needed for evaluating the helmothlz contributions
    temp_aux = eos.temperature_aux(T)
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]
    rhom = rho * Na

    dxhi00_drho = eos.dxhi00_drho
    diag_index = eos.diag_index
    xmi = x * eos.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    xhi00 = dxhi00_drho * rhom
    xhi = xhi00 * xm * np.matmul(xs, di03)
    
    # Monomer contribution
    xhi0, xhi1, xhi2, xhi3 = xhi
    xhi3_1 = 1-xhi3
    xhi23 = xhi2**3

    aHS = (xhi23/xhi3**2 - xhi0) * np.log(xhi3_1)
    aHS += 3 * xhi1 * xhi2 / xhi3_1
    aHS += xhi23 / xhi3 / xhi3_1**2
    aHS /= xhi0
    aMono = xm * aHS

    # Residual contribution
    aRes = eos.ares(x, rhom, T)

    tethai, Ai, Bi, Ci, Di = eos.viscosity_parameters.T
    aAux = (aRes - aMono)
    auxi = eos.ms/xm
    auxi *= Tad**tethai


    lnvisc_mix = np.dot(x, Ai)
    lnvisc_mix += np.dot(x * auxi, Bi) * aAux
    lnvisc_mix += np.dot(x * auxi**2, Ci) * aAux**2
    lnvisc_mix += np.dot(x * auxi**3, Di) * aAux**3

    visc = np.exp(lnvisc_mix) * visc_mixref * 1e3

    return visc
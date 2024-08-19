import numpy as np
from ..pcsaft_mixtures.ares import dares_dx
import types
from copy import deepcopy

kb = 1.3806488e-23  # [J/K] Boltzman's constant
Na = 6.022142e23  # [mol-1] Avogadro's Number
R = Na * kb    # [J mol-1 K-1] Ideal gas constant

def NETGP(saft, rho0NE = None, ksw = None):
    netsaft = deepcopy(saft)
    netsaft.logfugef_aux = types.MethodType(logfugef_NET, netsaft)
    netsaft.rho0NE = rho0NE
    netsaft.ksw = ksw
    return netsaft


def logfugef_NET(self, x, temp_aux, P, state, v0=None, Xass0=None):
    """
    logfugef_aux(x, temp_aux, P, state, v0, Xass0)

    Method that computes the chemical potential of the
    components in the mixture at given composition, temperature
    and pressure when in non-equilibrium.

    Parameters
    ----------
    x: array_like
        molar fraction array
    temp_aux : list
        temperature dependend parameters computed with temperature_aux(T)
    P: float
        pressure [Pa]
    state: string
        'L' for liquid phase and 'V' for vapour phase
    v0: float, optional
        initial guess for volume root [m^3/mol]
    Xass0: array, optional
        Initial guess for the calculation of fraction of non-bonded sites

    Returns
    -------
    mu: array_like
        chemical potential of the components
    v: float
        computed volume of the phase [m^3/mol]
    Xass : array
        computed fraction of nonbonded sites
    """

    beta = temp_aux[0]
    RT = Na/beta
    if v0 is None:
        rho0 = None
    else:
        rho0 = 1./v0

    rho, Xass = self.density_aux(x, temp_aux, P, state, rho0, Xass0)
    idx=1
    if state=="L":
        yi = self.yi
        idx_sol = yi!=0
        Pi = yi*P
        rho0NE = self.rho0p
        ksw = self.ksw
        xFeed = x.copy()
        xFeed[idx_sol] = 0
        xFeed /= xFeed.sum()
        rhopNE = NETPCSAFT(Pi, rho0NE, ksw, xFeed, nsw=2)     
        rho = rhopNE / (1 - x[idx_sol].sum())
        v = 1./rho
        Z = self.dP_drho_aux(x, rho, temp_aux)[0]/rho/RT
    else:
        
        v = 1./rho
        Z = P * v / RT
    rhomolecular = Na * rho
    
    ar, daresx, Xass = dares_dx(self, x, rhomolecular, temp_aux, Xass0)   
    mu = ar + (Z - 1.) + daresx - np.dot(x, daresx) + np.log(rho*RT)
    mu[idx] = -300 if state=="L" else mu[idx]
    #
    return mu, v, Xass
#we must supply the yi vector for the vapor phase
#for two solvents there must be some kind of feedback of the vapor phase compositions
#We give the composition of


def NETPCSAFT(Pi, rho0NE, ksw, xFeed, nsw=2):
    vFeed = xFeed / rho0NE / np.sum(xFeed/rho0NE)
    rhoF = np.sum(vFeed * rho0NE)
    kswmix = np.sum(vFeed * ksw)
    p_fun = 1 - np.sum(kswmix*(Pi)**nsw)
    rhoNE = rhoF * p_fun
    return rhoNE
from __future__ import division, print_function, absolute_import
import numpy as np
from  numba import njit

def association_config(eos):

    types = np.array(['B', 'P', 'N'])
    nozero = np.nonzero(eos.sites)
    types = types[nozero]
    ntypes = np.asarray(eos.sites)[nozero]
    nsites = len(types)
    S = np.array(eos.sites)
    S = S[S != 0]
    DIJ = np.zeros([nsites, nsites])
    int_i = []
    int_j = []
    for i in range(nsites):
        for j in range(nsites):
            bool1 = types[i] == 'B'
            bool2 = types[i] == 'P' and (types[j] == 'N' or types[j] == 'B')
            bool3 = types[i] == 'N' and (types[j] == 'P' or types[j] == 'B')
            if bool1 or bool2 or bool3:
                DIJ[i, j] = ntypes[j]
                int_i.append(i)
                int_j.append(j)

    indexabij = (int_i, int_j)
    diagasso = np.diag_indices(nsites)

    return S, DIJ, indexabij, nsites, diagasso

@njit
def Xass_solver(nsites, KIJ, diagasso, Xass):


    omega = 0.2

    for i in range(5):
        fo = 1. / (1. + KIJ@Xass)
        dXass = (1 - omega) * (fo - Xass)
        Xass += dXass

    for i in range(15):
        KIJXass = KIJ@Xass
        dQ = (1./Xass - 1.) - KIJXass
        HIJ = - 1. * KIJ
        HIJ -= np.diag((1. + KIJXass)/Xass)
        dXass = np.linalg.solve(HIJ, -dQ)
        Xass += dXass
        sucess = np.linalg.norm(dXass) < 1e-9
        if sucess:
            break

    return Xass

@njit
def dXass_drho(rho, Xass, DIJ, Dabij, dDabij_drho, CIJ):
    brho = -(DIJ*(Dabij + rho * dDabij_drho))@Xass
    brho *= Xass**2
    dXass = np.linalg.solve(CIJ, brho)
    return dXass

@njit
def d2Xass_drho(rho, Xass, dXass_drho, DIJ, Dabij, dDabij_drho, d2Dabij_drho,
                CIJ):

    b2rho = (2 * DIJ * (Dabij + rho * dDabij_drho))@dXass_drho
    b2rho += ( DIJ * (2 * dDabij_drho + rho * d2Dabij_drho))@Xass
    b2rho *= -Xass**2
    b2rho += (2*dXass_drho**2)/Xass

    d2Xass_drho = np.linalg.solve(CIJ, b2rho)
    return d2Xass_drho

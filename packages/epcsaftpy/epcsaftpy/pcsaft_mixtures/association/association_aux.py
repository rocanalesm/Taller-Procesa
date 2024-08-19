from __future__ import division, print_function, absolute_import
import numpy as np
from scipy.optimize import fsolve
from  numba import njit

def association_config(sitesmix, eos):

    compindex = []
    for i, j in enumerate(sitesmix):
        compindex += np.count_nonzero(j) * [i]
    compindex = np.asarray(compindex)

    types = np.array([['B', 'P', 'N']])
    n = len(sitesmix)
    n = eos.nc
    types = np.repeat(types, n, axis=0)
    nozero = np.nonzero(sitesmix)
    types = types[nozero]
    ntypes = np.asarray(sitesmix)[nozero]
    nsites = len(types)

    S = np.array(sitesmix)
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
    indexab = (compindex[int_i], compindex[int_j])
    diagasso = np.diag_indices(nsites)

    dxjdx = np.zeros([eos.nc, nsites])
    for i in range(eos.nc):
        dxjdx[i] = compindex == i

    return S, DIJ, compindex, indexabij, indexab, nsites, dxjdx, diagasso

@njit
def fobj_xass(Xass, xj, aux_asso, diagasso):
    fo = Xass - 1 / (1 + aux_asso@(xj*Xass))
    return fo

@njit
def fobj_xass_jac(Xass, xj, aux_asso, diagasso):
    den = 1 + aux_asso@(xj*Xass)
    dfo = ((aux_asso*xj).T/den**2).T
    dfo += np.eye(len(xj))
    return dfo

@njit
def fobj_xass_aux(Xass, xj, aux_asso, diagasso):
    den = 1 + aux_asso@(xj*Xass)
    fo = Xass - 1 / den
    dfo = ((aux_asso*xj).T/den**2).T
    dfo += np.eye(len(xj))
    return fo, dfo


@njit
def Xass_solver(nsites, xj, rho, DIJ, Dabij, diagasso, Xass0):
    
    dtype_cal = Xass0**0
    xj = xj * dtype_cal
    DIJ = DIJ * dtype_cal
    
    Xass = 1. * Xass0
    aux_asso = rho * DIJ * Dabij

    omega = 0.2
    for i in range(5):
        den = 1. + aux_asso@(xj*Xass)
        fo = 1. / den
        dXass = (1 - omega) * (fo - Xass)
        Xass += dXass


    KIJ = np.outer(xj, xj) * aux_asso

    KIJXass = KIJ@Xass
    dQ = xj * (1/Xass - 1) - KIJXass
    HIJ = -1 * KIJ
    HIJ -= np.diag((xj + KIJXass)/Xass)
    for i in range(15):
        dXass = np.linalg.solve(HIJ, -dQ)
        Xnew = Xass + dXass

        is_nan = np.isnan(Xnew.real)
        Xnew[is_nan] = 0.2

        Xnew_neg = Xnew.real < 0
        Xnew[Xnew_neg] = 0.2*Xass[Xnew_neg]
        Xass = Xnew
        KIJXass = KIJ@Xass
        dQ = xj * (1/Xass - 1) - KIJXass
        sucess = np.linalg.norm(dQ.real) < 1e-10
        if sucess:
            break
        HIJ = -1 * KIJ
        HIJ -= np.diag((xj + KIJXass)/Xass)


    return Xass

def Xass_solver_aux(nsites, xj, rho, DIJ, Dabij, diagasso, Xass0):
    Xass = 1. * Xass0
    aux_asso = rho * DIJ * Dabij

    omega = 0.2
    for i in range(5):
        den = 1. + aux_asso@(xj*Xass)
        fo = 1. / den
        dXass = (1 - omega) * (fo - Xass)
        Xass += dXass
        
    if not isinstance(Xass[0], np.complex128):
        Xass = fsolve(fobj_xass, x0=Xass, args=(xj, aux_asso, diagasso),
                      fprime=fobj_xass_jac)
    else:
        fo, dfo = fobj_xass_aux(Xass, xj, aux_asso, diagasso)
        for i in range(15):
            dX = np.matmul(np.linalg.inv(dfo),fo)
            Xass -= dX
            if np.linalg.norm(dX.real) < 1e-6:
                break
            fo, dfo = fobj_xass_aux(Xass, xj, aux_asso, diagasso)

    return Xass

"""
def Xass_solver(nsites, xj, rho, DIJ, Dabij, diagasso, Xass0):

    Xass = 1. * Xass0
    aux_asso = rho * DIJ * Dabij

    omega = 0.2
    for i in range(5):
        den = 1. + aux_asso@(xj*Xass)
        fo = 1. / den
        dXass = (1 - omega) * (fo - Xass)
        Xass += dXass

    bool_method = not np.any(xj == 0.)

    if bool_method:
        KIJ = np.outer(xj, xj) * aux_asso

        KIJXass = KIJ@Xass
        dQ = xj * (1/Xass - 1) - KIJXass
        HIJ = -1 * KIJ
        HIJ[diagasso] -= (xj + KIJXass)/Xass
        for i in range(15):
            dXass = np.linalg.solve(HIJ, -dQ)
            Xnew = Xass + dXass

            is_nan = np.isnan(Xnew)
            Xnew[is_nan] = 0.2

            Xnew_neg = Xnew < 0
            Xnew[Xnew_neg] = 0.2*Xass[Xnew_neg]
            Xass = Xnew
            KIJXass = KIJ@Xass
            dQ = xj * (1/Xass - 1) - KIJXass
            sucess = np.linalg.norm(dQ) < 1e-10
            if sucess:
                break
            HIJ = -1 * KIJ
            HIJ[diagasso] -= (xj + KIJXass)/Xass
    else:
        Xass = fsolve(fobj_xass, x0=Xass, args=(xj, aux_asso, diagasso),
                      fprime=fobj_xass_jac)

    return Xass
"""

@njit
def CIJ_matrix(rhom, xj, Xass, DIJ, Dabij, diagasso):
    CIJ = rhom * np.outer(Xass**2, xj) * Dabij * DIJ
    CIJ += np.eye(len(xj))
    return CIJ

@njit
def dXass_drho(rhom, xj, Xass, DIJ, Dabij, dDabij_drho, CIJ):
    brho = -(DIJ*(Dabij + rhom * dDabij_drho))@(xj * Xass)
    brho *= Xass**2
    dXass = np.linalg.solve(CIJ, brho)
    return dXass

@njit
def d2Xass_drho(rhom, xj, Xass, dXass_drho, DIJ, Dabij, dDabij_drho,
                d2Dabij_drho, CIJ):

    b2rho = (2 * DIJ * (Dabij + rhom * dDabij_drho))@(xj * dXass_drho)
    b2rho += ( DIJ * (2 * dDabij_drho + rhom * d2Dabij_drho))@(xj * Xass)
    b2rho *= -Xass**2
    b2rho += (2*dXass_drho**2)/Xass


    d2Xass = np.linalg.solve(CIJ, b2rho)
    return d2Xass

def dXass_dx(rhom, xj, Xass, DIJ, Dabij, dDabij_dx, dxjdx, CIJ):

    bx = (dxjdx * Xass)@(DIJ*Dabij).T + (DIJ * dDabij_dx)@(xj * Xass)
    bx *= Xass**2
    bx *= -rhom
    dXass = np.linalg.solve(CIJ, bx.T)

    return dXass.T

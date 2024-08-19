from __future__ import division, print_function, absolute_import
import numpy as np
from ..constants import Na

def dPsaft_fun(rho, temp_aux, saft):
    rhomolecular = Na * rho
    global Xass
    da, Xass = saft.d2afcn_aux(rhomolecular, temp_aux, Xass)
    afcn, dafcn, d2afcn = da
    dPsaft = 2 * rhomolecular * dafcn + rhomolecular**2 * d2afcn
    return dPsaft


def Psaft_obj(rho, temp_aux, saft, Pspec):
    rhomolecular = Na * rho
    global Xass
    da, Xass = saft.dafcn_aux(rhomolecular, temp_aux, Xass)
    afcn, dafcn = da
    Psaft = rhomolecular**2 * dafcn / Na
    return Psaft - Pspec


def density_newton(rho0, temp_aux, P, Xass, rhomax, state, saft):
    rho = 1.*rho0
    for i in range(15):
        Psaft, dPsaft, Xass = saft.dP_drho_aux(rho, temp_aux, Xass)
        while dPsaft<=0:
            if state == 'L':
                rho = (rhomax + rho)/2
            elif state == 'V':
                rho = rho/2
            Psaft, dPsaft, Xass = saft.dP_drho_aux(rho, temp_aux, Xass)

        FO = Psaft - P
        dFO = dPsaft
        drho = FO/dFO
        rho_new = rho - drho
        if rho_new <= 0:
            rho = rho/2
        elif rho_new > rhomax:
            rho = (rhomax + rho)/2
        else:
            rho = np.copy(rho_new)
        if np.abs(drho/rho) < 1e-12:
            break
    return rho, Xass

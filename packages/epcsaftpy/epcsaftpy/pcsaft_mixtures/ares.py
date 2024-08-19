from __future__ import division, print_function, absolute_import
import numpy as np

from .monomer.aHS import ahs, dahs_dxhi00, d2ahs_dxhi00,  dahs_dx, dahs_dxxhi
from .chain.gHS import ghs, dghs_dxhi00, d2ghs_dxhi00,  dghs_dx, dghs_dxxhi
from .chain.achain import achain, dachain_dxhi00, d2achain_dxhi00,  dachain_dx, dachain_dxxhi
from .dispersion.adisp import adisp, dadisp_dxhi00, d2adisp_dxhi00, dadisp_dx, dadisp_dxxhi

from .association.association_aux import Xass_solver, Xass_solver_aux, CIJ_matrix
from .association.association_aux import dXass_drho, d2Xass_drho, dXass_dx 

from .polar.apolarJC import apolarJC, dapolarJC_drho, d2apolarJC_drho
from .polar.apolarJC import d2apolarJC_dx, d2apolarJC_dxrho

from .debye_huckel.aDH import adh, dadh_drho, d2adh_drho, dadh_dx, dadh_dxrho
from .born.aBorn import aborn, daborn_drho, d2aborn_drho, daborn_dx, daborn_dxrho


def xhi_eval(xhi00, xs, xmi, xm, di03):
    xhi = xhi00 * xm * np.matmul(xs, di03)
    dxhi_dxhi00 = np.matmul(xmi, di03)
    dxhi_dxhi00[0] = xm
    return xhi, dxhi_dxhi00


def dxhi_dx_eval(xhi00, xs, xmi, xm, ms, di03):
    xhi = xhi00 * xm * np.matmul(xs, di03)
    dxhi_dxhi00 = np.matmul(xmi, di03)
    dxhi_dxhi00[0] = xm
    dxhi_dx = (xhi00 * di03.T * ms)
    return xhi, dxhi_dxhi00, dxhi_dx

def mix_er(x, Mw_solv, index_ions, index_solv, eri_ions, eri_solv):
    x_ions = x[index_ions]
    x_solv = x[index_solv]
    Mwx = x_solv * Mw_solv
    sum_Mwx = np.sum(Mwx)
    w_solv = Mwx / sum_Mwx
    
    er_ions = np.dot(eri_ions, x_ions)
    er_solv = np.dot(eri_solv, w_solv)
    er = er_ions + np.dot(er_solv, np.sum(x_solv))

    return er

def mix_derx(x, Mw_solv, index_ions, index_solv, eri_ions, eri_solv):
    x_ions = x[index_ions]
    x_solv = x[index_solv]
    Mwx = x_solv * Mw_solv
    sum_Mwx = np.sum(Mwx)
    w_solv = Mwx / sum_Mwx
    
    er_ions = np.dot(eri_ions, x_ions)
    er_solv = np.dot(eri_solv, w_solv)
    er = er_ions + np.dot(er_solv, np.sum(x_solv))
    
    dwdx = Mw_solv * np.eye(len(x_solv)) * sum_Mwx
    dwdx -= np.outer(Mwx, Mw_solv)
    dwdx /= sum_Mwx**2
    
    der_solv = eri_solv@dwdx
    derx = np.zeros_like(x)
    derx[index_ions] = eri_ions
    derx[index_solv] = er_solv + der_solv * np.sum(x_solv)
    return er, derx




def ares(self, x, rho, temp_aux, Xass0=None):
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]

    # Parameters needed for evaluating the helmothlz contributions
    dxhi00_drho = self.dxhi00_drho
    diag_index = self.diag_index
    xmi = x * self.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    xhi00 = dxhi00_drho * rho
    xhi, dxhi_dxhi00 = xhi_eval(xhi00, xs, xmi, xm, di03)
    
    # Monomer contribution
    aHS = ahs(xhi)
    aMono = xm * aHS

    # Chain contribution
    gHS = ghs(xhi, aux_dii, aux_dii2)
    aChain = achain(x, self.ms, gHS, diag_index)
    
    # Dispersion contribution
    eta = xhi[-1]
    aDips = adisp(xhi00, xm, xmi, eta, es3ij, e2s3ij)
    
    ares = aMono + aChain + aDips
    
    # Association contribution
    if self.assoc_bool:
        xj = x[self.compindex]
        
        Dab = aux_Dab * gHS
        if type(aHS) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
        else:
            Dabij = np.zeros([self.nsites, self.nsites])
        Dabij[self.indexabij] = Dab[self.indexab]
        
        if Xass0 is None:
            KIJ = rho * (self.DIJ*Dabij)
            sumD = np.sum(xj*KIJ, axis = 1)
            sumD[sumD==0] = 1e-12
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)
            is_nan = np.isnan(Xass0.real)
            Xass0[is_nan] = 0.2
        bool_method = not np.any(xj == 0.)
        if bool_method:
            Xass = Xass_solver(self.nsites, xj, rho, self.DIJ, Dabij,
                               self.diagasso, Xass0)
        else:
            Xass = Xass_solver_aux(self.nsites, xj, rho, self.DIJ, Dabij,
                                   self.diagasso, Xass0)
        ares += np.dot(self.S * xj, (np.log(Xass) - Xass/2 + 1/2))

    else:
        Xass = None

    # Polar contribution (JC)
    if self.polar_bool:
        aux1_polarij, aux2_polarijk, dij3, dijk3 = temp_aux[9:13]
        apolar = apolarJC(x, rho, eta, aux1_polarij, 
                          aux2_polarijk, dij3, dijk3)
        ares += apolar
        
    # Electrolite contribution
    if self.elec_bool:
        dii, betae2_e0 = temp_aux[13:15]
        
        # Mixing rules for er
        er = mix_er(x, self.Mw_solv, self.index_ions, self.index_solv, 
                    self.eri_ions, self.eri_solv)
        
        xz2 = x * self.z**2
        # Debye-Huckel contribution
        aDH = adh(rho, xz2, er, dii, betae2_e0)
        ares += aDH
        # Born contribution 
        aBorn = aborn(xz2, er, dii, betae2_e0)
        ares += aBorn  

    return ares, Xass


def dares_drho(self, x, rho, temp_aux, Xass0=None):
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]

    # Parameters needed for evaluating the helmothlz contributions
    dxhi00_drho = self.dxhi00_drho
    diag_index = self.diag_index
    xmi = x * self.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    xhi00 = dxhi00_drho * rho
    xhi, dxhi_dxhi00 = xhi_eval(xhi00, xs, xmi, xm, di03)

    # Monomer contribution
    aHS = dahs_dxhi00(xhi, dxhi_dxhi00)
    aMono = xm * aHS

    # Chain contribution
    gHS = dghs_dxhi00(xhi, dxhi_dxhi00, aux_dii, aux_dii2)
    aChain = dachain_dxhi00(x, self.ms, gHS, diag_index)
    
    # Dispersion contribution 
    eta = xhi[-1]
    deta_dxhi00 = dxhi_dxhi00[-1]
    aDisp = dadisp_dxhi00(xhi00, xm, xmi, eta, deta_dxhi00, es3ij, e2s3ij)


    ares = aMono + aChain + aDisp
    ares *= np.array([1, dxhi00_drho])
    
    # Association contribution
    if self.assoc_bool:
        xj = x[self.compindex]
        ghs_eval = gHS[0]
        dghs_rho = gHS[1]*dxhi00_drho
        
        Dab = aux_Dab * ghs_eval
        dDab_drho = Dab * dghs_rho/ghs_eval
        
        if type(aHS[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])
            
        Dabij[self.indexabij] = Dab[self.indexab]
        dDabij_drho[self.indexabij] = dDab_drho[self.indexab]
        
        if Xass0 is None:
            KIJ = rho * (self.DIJ*Dabij)
            sumD = np.sum(xj*KIJ, axis = 1)
            sumD[sumD==0] = 1e-12
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)
            is_nan = np.isnan(Xass0.real)
            Xass0[is_nan] = 0.2
            
        bool_method = not np.any(xj == 0.)
        if bool_method:
            Xass = Xass_solver(self.nsites, xj, rho, self.DIJ, Dabij,
                               self.diagasso, Xass0)
        else:
            Xass = Xass_solver_aux(self.nsites, xj, rho, self.DIJ, Dabij,
                                   self.diagasso, Xass0)
        CIJ = CIJ_matrix(rho, xj, Xass, self.DIJ, Dabij, self.diagasso)
        dXass = dXass_drho(rho, xj, Xass, self.DIJ, Dabij, dDabij_drho, CIJ)
        ares[0] += np.dot(self.S * xj, (np.log(Xass) - Xass/2 + 1/2))
        ares[1] += np.dot(self.S*xj, (1/Xass - 1/2) * dXass)

    else:
        Xass = None

    # Polar contribution (JC)
    if self.polar_bool:
        aux1_polarij, aux2_polarijk, dij3, dijk3 = temp_aux[9:13]
        dapolar = dapolarJC_drho(x, rho, eta, aux1_polarij, aux2_polarijk, 
                                  dij3, dijk3)
        ares += dapolar
        
    # Electrolite contribution
    if self.elec_bool:
        dii, betae2_e0 = temp_aux[13:15]
        
        # Mixing rule for er
        er = mix_er(x, self.Mw_solv, self.index_ions, self.index_solv, 
                    self.eri_ions, self.eri_solv)
        
        xz2 = x * self.z**2
        # Debye-Huckel contribution
        aDH = dadh_drho(rho, xz2, er, dii, betae2_e0)
        ares += aDH
        # Born contribution 
        aBorn = daborn_drho(xz2, er, dii, betae2_e0)
        ares += aBorn

    return ares, Xass


def d2ares_drho(self, x, rho, temp_aux, Xass0=None):
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]

    # Parameters needed for evaluating the helmothlz contributions
    dxhi00_drho = self.dxhi00_drho
    diag_index = self.diag_index
    xmi = x * self.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    xhi00 = dxhi00_drho * rho

    xhi, dxhi_dxhi00 = xhi_eval(xhi00, xs, xmi, xm, di03)

    # Monomer contribution
    aHS = d2ahs_dxhi00(xhi, dxhi_dxhi00)
    aMono = xm * aHS

    # Chain contribution
    gHS = d2ghs_dxhi00(xhi, dxhi_dxhi00, aux_dii, aux_dii2)
    aChain = d2achain_dxhi00(x, self.ms, gHS, diag_index)
    
    # Dispersion contribution
    eta = xhi[-1]
    deta_dxhi00 = dxhi_dxhi00[-1]
    aDisp = d2adisp_dxhi00(xhi00, xm, xmi, eta, deta_dxhi00, es3ij, e2s3ij)
    
    ares = aMono + aChain + aDisp
    ares *= np.array([1., dxhi00_drho, dxhi00_drho**2])
    
    # Association contribution
    if self.assoc_bool:
        xj = x[self.compindex]
            
        ghs_eval = gHS[0]
        dghs_rho = gHS[1]*dxhi00_drho
        d2ghs_rho = gHS[2]*dxhi00_drho**2
        
        Dab = aux_Dab * ghs_eval
        dDab_drho = Dab * dghs_rho/ghs_eval
        d2Dab_drho = Dab * d2ghs_rho/ghs_eval
        
        if type(aHS[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            d2Dabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')

        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])
            d2Dabij_drho = np.zeros([self.nsites, self.nsites])

        Dabij[self.indexabij] = Dab[self.indexab]
        dDabij_drho[self.indexabij] = dDab_drho[self.indexab]
        d2Dabij_drho[self.indexabij] = d2Dab_drho[self.indexab]
        
        if Xass0 is None:
            KIJ = rho * (self.DIJ*Dabij)
            sumD = np.sum(xj*KIJ, axis = 1)
            sumD[sumD==0] = 1e-12
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)
            is_nan = np.isnan(Xass0.real)
            Xass0[is_nan] = 0.2
            
        bool_method = not np.any(xj == 0.)
        if bool_method:
            Xass = Xass_solver(self.nsites, xj, rho, self.DIJ, Dabij,
                               self.diagasso, Xass0)
        else:
            Xass = Xass_solver_aux(self.nsites, xj, rho, self.DIJ, Dabij,
                                   self.diagasso, Xass0)
        CIJ = CIJ_matrix(rho, xj, Xass, self.DIJ, Dabij, self.diagasso)
        dXass = dXass_drho(rho, xj, Xass, self.DIJ, Dabij, dDabij_drho,
                           CIJ)
        d2Xass = d2Xass_drho(rho, xj, Xass, dXass, self.DIJ, Dabij,
                             dDabij_drho, d2Dabij_drho, CIJ)

        aux1 = np.log(Xass) - Xass/2 + 1/2
        aux2 = 1/Xass - 1/2

        ares[0] += np.dot(self.S*xj, aux1)
        ares[1] += np.dot(self.S*xj, aux2 * dXass)
        ares[2] += np.dot(self.S*xj, -(dXass/Xass)**2+d2Xass*aux2)

    else:
        Xass = None
        
    # Polar contribution (JC)
    if self.polar_bool:
        aux1_polarij, aux2_polarijk, dij3, dijk3 = temp_aux[9:13]
        dapolar = d2apolarJC_drho(x, rho, eta, aux1_polarij, aux2_polarijk, 
                                  dij3, dijk3)
        ares += dapolar
        
    # Electrolite contribution
    if self.elec_bool:
        dii, betae2_e0 = temp_aux[13:15]
        
        # Mixing rule for er
        er = mix_er(x, self.Mw_solv, self.index_ions, self.index_solv, 
                    self.eri_ions, self.eri_solv)
        
        xz2 = x * self.z**2
        # Debye-Huckel contribution
        aDH = d2adh_drho(rho, xz2, er, dii, betae2_e0)
        ares += aDH
        # Born contribution 
        aBorn = d2aborn_drho(xz2, er, dii, betae2_e0)
        ares += aBorn

    return ares, Xass


def dares_dx(self, x, rho, temp_aux, Xass0=None):
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]

    # Parameters needed for evaluating the helmothlz contributions    
    dxhi00_drho = self.dxhi00_drho
    diag_index = self.diag_index

    xmi = x * self.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    dxs_dx = - np.multiply.outer(self.ms, self.ms * x) / xm**2
    dxs_dx[diag_index] += self.ms / xm
    xhi00 = dxhi00_drho * rho
    out = dxhi_dx_eval(xhi00, xs, xmi, xm, self.ms, di03)
    xhi, dxhi_dxhi00, dxhi_dx = out

  
    # Monomer contribution
    aHS, daHSx = dahs_dx(xhi, dxhi_dx)
    aMono = xm * aHS
    daMonox = self.ms * aHS + xm * daHSx

    # Chain contribution
    gHS, dgHSx = dghs_dx(xhi, dxhi_dx, aux_dii, aux_dii2)
    aChain, daChainx = dachain_dx(x, self.ms, gHS, dgHSx, diag_index)
    
    # Dispersion contribution
    eta = xhi[-1]
    deta_dx = dxhi_dx[-1]
    aDisp, daDispx = dadisp_dx(x, xhi00, xm, xmi, self.ms, eta, 
                               deta_dx, es3ij, e2s3ij, mes3ij, me2s3ij)
    
    ares = aMono + aChain + aDisp
    daresx = daMonox + daChainx + daDispx

    # Association contribution
    if self.assoc_bool:
        xj = x[self.compindex]
        
        Dab = aux_Dab * gHS
        dDab_dx = aux_Dab * dgHSx
    
        Dabij = np.zeros([self.nsites, self.nsites])
        dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites])
        
        if type(aHS) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites], dtype = 'complex_')

        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites])
        
        Dabij[self.indexabij] = Dab[self.indexab]
        dDabij_dx[:, self.indexabij[0], self.indexabij[1]] = dDab_dx[:, self.indexab[0], self.indexab[1]]
        
        if Xass0 is None:
            KIJ = rho * (self.DIJ*Dabij)
            sumD = np.sum(xj*KIJ, axis = 1)
            sumD[sumD==0] = 1e-12
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)
            is_nan = np.isnan(Xass0.real)
            Xass0[is_nan] = 0.2
            
        bool_method = not np.any(xj == 0.)
        if bool_method:
            Xass = Xass_solver(self.nsites, xj, rho, self.DIJ, Dabij,
                               self.diagasso, Xass0)
        else:
            Xass = Xass_solver_aux(self.nsites, xj, rho, self.DIJ, Dabij,
                                   self.diagasso, Xass0)
        CIJ = CIJ_matrix(rho, xj, Xass, self.DIJ, Dabij, self.diagasso)
        dXassx = dXass_dx(rho, xj, Xass, self.DIJ, Dabij, dDabij_dx,
                          self.dxjdx, CIJ)

        aux1 = np.log(Xass) - Xass/2 + 1/2
        aux2 = 1/Xass - 1/2

        aasso = np.dot(self.S*xj, aux1)
        daassox = (self.dxjdx * aux1 + dXassx * xj * aux2)@self.S
        ares += aasso
        daresx += daassox


    else:
        Xass = None

    # Polar contribution (JC)
    if self.polar_bool:
        aux1_polarij, aux2_polarijk, dij3, dijk3 = temp_aux[9:13]
        a, dax = d2apolarJC_dx(x, rho, eta, aux1_polarij, aux2_polarijk, 
                               dij3, dijk3, deta_dx)
        ares += a
        daresx += dax
        
    # Electrolite contribution
    if self.elec_bool:
        dii, betae2_e0 = temp_aux[13:15]
        
        # Mixing rules for er
        er, derx = mix_derx(x, self.Mw_solv, self.index_ions, self.index_solv, 
                            self.eri_ions, self.eri_solv)

        xz2 = x * self.z**2
        # Debye-Huckel contribution
        aDH, daDHx = dadh_dx(rho, derx, self.z, xz2, er, dii, betae2_e0)
        ares += aDH
        daresx += daDHx
        # Born contribution 
        aBorn, daBornx = daborn_dx(xz2, self.z, er, derx, dii, betae2_e0)
        ares += aBorn
        daresx += daBornx

    return ares, daresx, Xass


def dares_dxrho(self, x, rho, temp_aux, Xass0=None):
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]

    # Parameters needed for evaluating the helmothlz contributions    

    dxhi00_drho = self.dxhi00_drho
    diag_index = self.diag_index
    xmi = x * self.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    dxs_dx = - np.multiply.outer(self.ms, self.ms * x) / xm**2
    dxs_dx[diag_index] += self.ms / xm
    xhi00 = dxhi00_drho * rho
    xhi, dxhi_dxhi00, dxhi_dx = dxhi_dx_eval(xhi00, xs, xmi, xm, self.ms, di03)


    # Monomer contribution
    aHS, daHSx = dahs_dxxhi(xhi, dxhi_dxhi00, dxhi_dx)
    aMono = xm * aHS
    daMonox = self.ms * aHS[0] + xm * daHSx

    # Chain contribution
    gHS, dgHSx = dghs_dxxhi(xhi, dxhi_dxhi00, dxhi_dx, aux_dii, aux_dii2)
    aChain, daChainx = dachain_dxxhi(x, self.ms, gHS, dgHSx, diag_index)
    
    # Dispersion contribution
    eta = xhi[-1]
    deta_dxhi00 = dxhi_dxhi00[-1]
    deta_dx = dxhi_dx[-1]
    aDisp, daDispx = dadisp_dxxhi(x, xhi00, xm, xmi, self.ms, eta, deta_dx, 
                     deta_dxhi00, es3ij, e2s3ij, mes3ij, me2s3ij)

    ares = aMono + aChain + aDisp
    daresx = daMonox + daChainx + daDispx

    ares *= np.array([1, dxhi00_drho])

    # Association contribution
    if self.assoc_bool:
        xj = x[self.compindex]
        
        ghs_eval = gHS[0]
        dghs_rho = gHS[1]*dxhi00_drho
        
        Dab = aux_Dab * ghs_eval
        dDab_drho = aux_Dab * dghs_rho
        dDab_dx = aux_Dab * dgHSx
                
        if type(aHS[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites], dtype = 'complex_')

        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])
            dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites])

        Dabij[self.indexabij] = Dab[self.indexab]
        dDabij_drho[self.indexabij] = dDab_drho[self.indexab]
        dDabij_dx[:, self.indexabij[0], self.indexabij[1]] = dDab_dx[:, self.indexab[0], self.indexab[1]]
        
        if Xass0 is None:
            KIJ = rho * (self.DIJ*Dabij)
            sumD = np.sum(xj*KIJ, axis = 1)
            sumD[sumD==0] = 1e-12
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)
            is_nan = np.isnan(Xass0.real)
            Xass0[is_nan] = 0.2
            
        bool_method = not np.any(xj == 0.)
        if bool_method:
            Xass = Xass_solver(self.nsites, xj, rho, self.DIJ, Dabij,
                               self.diagasso, Xass0)
        else:
            Xass = Xass_solver_aux(self.nsites, xj, rho, self.DIJ, Dabij,
                                   self.diagasso, Xass0)
        CIJ = CIJ_matrix(rho, xj, Xass, self.DIJ, Dabij, self.diagasso)
        dXassx = dXass_dx(rho, xj, Xass, self.DIJ, Dabij, dDabij_dx,
                          self.dxjdx, CIJ)
        dXass = dXass_drho(rho, xj, Xass, self.DIJ, Dabij, dDabij_drho, CIJ)

        aux1 = np.log(Xass) - Xass/2 + 1/2
        aux2 = 1/Xass - 1/2

        aasso = np.dot(self.S*xj, aux1)
        daasso = np.dot(self.S*xj, aux2 * dXass)
        

        ares[0] += aasso
        ares[1] += daasso


        daassox = (self.dxjdx * aux1 + dXassx * xj * aux2)@self.S
        daresx += daassox
    else:
        Xass = None
        
    # Polar contribution (JC)
    if self.polar_bool:
        aux1_polarij, aux2_polarijk, dij3, dijk3 = temp_aux[9:13]
        a, dax = d2apolarJC_dxrho(x, rho, eta, aux1_polarij, aux2_polarijk, 
                                  dij3, dijk3, deta_dx)
        ares += a
        daresx += dax
        
    # Electrolite contribution
    if self.elec_bool:
        dii, betae2_e0 = temp_aux[13:15]
        
        # Mixing rule for er
        er, derx = mix_derx(x, self.Mw_solv, self.index_ions, self.index_solv, 
                            self.eri_ions, self.eri_solv)

        xz2 = x * self.z**2
        # Debye-Huckel contribution
        aDH, daDHx = dadh_dxrho(rho, derx, self.z, xz2, er, dii, betae2_e0)
        ares += aDH
        daresx += daDHx
        # Born contribution 
        aBorn, daBornx = daborn_dxrho(xz2, self.z, er, derx, dii, betae2_e0)
        ares += aBorn
        daresx += daBornx

    return ares, daresx, Xass


def all_contributions(self, x, rho, temp_aux, Xass0=None):
    beta, di03, aux_dii, aux_dii2 = temp_aux[0:4]
    es3ij, e2s3ij, mes3ij, me2s3ij, aux_Dab = temp_aux[4:9]

    # Parameters needed for evaluating the helmothlz contributions    
    dxhi00_drho = self.dxhi00_drho
    diag_index = self.diag_index
    xmi = x * self.ms
    xm = np.sum(xmi)
    xs = xmi / xm
    dxs_dx = - np.multiply.outer(self.ms, self.ms * x) / xm**2
    dxs_dx[diag_index] += self.ms / xm
    xhi00 = dxhi00_drho * rho
    xhi, dxhi_dxhi00, dxhi_dx = dxhi_dx_eval(xhi00, xs, xmi, xm, self.ms, di03)


    # Monomer contribution
    aHS, daHSx = dahs_dxxhi(xhi, dxhi_dxhi00, dxhi_dx)
    aHS = d2ahs_dxhi00(xhi, dxhi_dxhi00)
    aMono = xm * aHS
    daMonox = self.ms * aHS[0] + xm * daHSx

    # Chain contribution
    gHS, dgHSx = dghs_dxxhi(xhi, dxhi_dxhi00, dxhi_dx, aux_dii, aux_dii2)
    aChain, daChainx = dachain_dxxhi(x, self.ms, gHS, dgHSx, diag_index)
    gHS = d2ghs_dxhi00(xhi, dxhi_dxhi00, aux_dii, aux_dii2)
    aChain = d2achain_dxhi00(x, self.ms, gHS, diag_index)
    
    
    # Dispersion contribution
    eta = xhi[-1]
    deta_dxhi00 = dxhi_dxhi00[-1]
    deta_dx = dxhi_dx[-1]
    aDisp, daDispx = dadisp_dxxhi(x, xhi00, xm, xmi, self.ms, eta, deta_dx, 
                     deta_dxhi00, es3ij, e2s3ij, mes3ij, me2s3ij)
    aDisp = d2adisp_dxhi00(xhi00, xm, xmi, eta, deta_dxhi00, es3ij, e2s3ij)

    aMono *= np.array([1., dxhi00_drho, dxhi00_drho**2])
    aChain *= np.array([1., dxhi00_drho, dxhi00_drho**2])
    aDisp *= np.array([1., dxhi00_drho, dxhi00_drho**2])

    # Association contribution
    aAsso = np.zeros_like(aDisp)
    daAssox = np.zeros_like(daDispx)
    if self.assoc_bool:
        xj = x[self.compindex]

        ghs_eval = gHS[0] * 1.
        dghs_rho = gHS[1] * dxhi00_drho
        d2ghs_rho = gHS[2] * dxhi00_drho**2
        

        
        Dab = aux_Dab * ghs_eval
        dDab_drho = aux_Dab * dghs_rho
        d2Dab_drho = aux_Dab * d2ghs_rho
        dDab_dx = aux_Dab * dgHSx
                
        if type(aHS[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            d2Dabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites], dtype = 'complex_')
        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])
            d2Dabij_drho = np.zeros([self.nsites, self.nsites])
            dDabij_dx = np.zeros([self.nc, self.nsites, self.nsites])

        Dabij[self.indexabij] = Dab[self.indexab]
        dDabij_drho[self.indexabij] = dDab_drho[self.indexab]
        d2Dabij_drho[self.indexabij] = d2Dab_drho[self.indexab]
        dDabij_dx[:, self.indexabij[0], self.indexabij[1]] = dDab_dx[:, self.indexab[0], self.indexab[1]]
        
        if Xass0 is None:
            KIJ = rho * (self.DIJ*Dabij)
            sumD = np.sum(xj*KIJ, axis = 1)
            sumD[sumD==0] = 1e-12
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)
            is_nan = np.isnan(Xass0.real)
            Xass0[is_nan] = 0.2
            
        bool_method = not np.any(xj == 0.)
        if bool_method:
            Xass = Xass_solver(self.nsites, xj, rho, self.DIJ, Dabij,
                               self.diagasso, Xass0)
        else:
            Xass = Xass_solver_aux(self.nsites, xj, rho, self.DIJ, Dabij,
                                   self.diagasso, Xass0)
        
        CIJ = CIJ_matrix(rho, xj, Xass, self.DIJ, Dabij, self.diagasso)
        dXass = dXass_drho(rho, xj, Xass, self.DIJ, Dabij, dDabij_drho,
                           CIJ)
        d2Xass = d2Xass_drho(rho, xj, Xass, dXass, self.DIJ, Dabij,
                             dDabij_drho, d2Dabij_drho, CIJ)
        dXassx = dXass_dx(rho, xj, Xass, self.DIJ, Dabij, dDabij_dx,
                          self.dxjdx, CIJ)
        


        aux1 = np.log(Xass) - Xass/2 + 1/2
        aux2 = 1/Xass - 1/2
        aa1 = np.dot(self.S*xj, aux1)
        aa2 = np.dot(self.S*xj, aux2 * dXass)
        aa3 = np.dot(self.S*xj, -(dXass/Xass)**2 + d2Xass*aux2)
        aAsso = np.array([aa1, aa2, aa3])
        
        daAssox = (self.dxjdx * aux1 + dXassx * xj * aux2)@self.S

        
        
    else:
        Xass = None
        
    # Polar contribution (JC)
    aPolar = np.zeros_like(aDisp)
    daPolarx = np.zeros_like(daDispx)
    if self.polar_bool:
        aux1_polarij, aux2_polarijk, dij3, dijk3 = temp_aux[9:13]
        _, daPolarx = d2apolarJC_dx(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3, deta_dx)
        aPolar = d2apolarJC_drho(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3)
    
    # Electrolite contribution
    aDH = np.zeros_like(aDisp)
    daDHx = np.zeros_like(daDispx)
    aBorn = np.zeros_like(aDisp)
    daBornx = np.zeros_like(daDispx)
    if self.elec_bool:
        dii, betae2_e0 = temp_aux[13:15]
        
        # Mixing rule for er
        er, derx = mix_derx(x, self.Mw_solv, self.index_ions, self.index_solv, 
                            self.eri_ions, self.eri_solv)

        xz2 = x * self.z**2
        # Debye-Huckel contribution
        _, daDHx = dadh_dxrho(rho, derx, self.z, xz2, er, dii, betae2_e0)
        aDH = d2adh_drho(rho, xz2, er, dii, betae2_e0)

        # Born contribution 
        _, daBornx = daborn_dxrho(xz2, self.z, er, derx, dii, betae2_e0)
        aBorn = d2aborn_drho(xz2, er, dii, betae2_e0)

        
    Mono = aMono, daMonox
    Chain = aChain, daChainx
    Disp = aDisp, daDispx
    Asso = aAsso, daAssox
    Polar = aPolar, daPolarx
    DH = aDH, daDHx
    Born = aBorn, daBornx
    return Mono, Chain, Disp,  Asso, Polar, DH, Born

from __future__ import division, print_function, absolute_import
import numpy as np


from .monomer.aHS import ahs, dahs_deta, d2ahs_deta
from .chain.achain import achain, dachain_deta, d2achain_deta
from .chain.gHS import ghs, dghs_deta, d2ghs_deta
from .dispersion.adisp import adisp, dadisp_deta, d2adisp_deta
from .polar.apolarJC import apolarJC, dapolarJC_deta, d2apolarJC_deta

from .association.association_aux import Xass_solver, dXass_drho, d2Xass_drho



def ares(self, rho, temp_aux, Xass0=None):
    # Parameters needed for evaluating the helmothlz contributions
    beta, eps_beta, d, dia3, mes3, m2e2s3, deta_drho, Fab, betaxpmu2 = temp_aux
    eta = deta_drho * rho
    
    # Monomer contribution
    a_hs = ahs(eta)

    # Chain contribution
    ghs_eval = ghs(eta)
    a_chain = achain(self.ms, ghs_eval)

    # Dispersion contribution
    a_disp = adisp(self.ms, eta, dia3, self.ai, self.bi, m2e2s3, mes3)
    
    # Polar contribution
    apolar = 0
    if self.polar_bool:
        apolar = apolarJC(eta, self.ms, dia3, betaxpmu2)

    # Total Helmolthz
    a = (self.ms*a_hs + a_chain) + a_disp + apolar

    if self.assoc_bool:
        Dab = self.sigma3 * Fab * self.kappaABij * ghs_eval
        if type(a_hs) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
        else:
            Dabij = np.zeros([self.nsites, self.nsites])
        Dabij[self.indexabij] = Dab
        
        KIJ = rho * (self.DIJ*Dabij)
        if Xass0 is None:
            sumD = np.sum(KIJ, axis = 1)
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)

        Xass = Xass_solver(self.nsites, KIJ, self.diagasso, Xass0)
        a += np.dot(self.S, (np.log(Xass) - Xass/2 + 1/2))
    else:
        Xass = None


    return a, Xass


def dares_drho(self, rho, temp_aux, Xass0=None):
    # Parameters needed for evaluating the helmothlz contributions
    beta, eps_beta, d, dia3, mes3, m2e2s3, deta_drho, Fab, betaxpmu2 = temp_aux
    eta = deta_drho * rho
    
    # Monomer contribution
    a_hs = dahs_deta(eta)

    # Chain contribution
    dghs = dghs_deta(eta)
    a_chain = dachain_deta(self.ms, dghs)

    # Dispersion contribution
    a_disp = dadisp_deta(self.ms, eta, dia3, self.ai, self.bi, m2e2s3, mes3)

    # Polar contribution
    apolar = np.zeros(2)
    if self.polar_bool:
        apolar = dapolarJC_deta(eta, self.ms, dia3, betaxpmu2)

    # Total Helmolthz
    a = (self.ms*a_hs + a_chain) + a_disp + apolar
    a *= np.array([1, deta_drho]) 


    if self.assoc_bool:

        ghs_eval = dghs[0]
        dghs_rho = dghs[1]*deta_drho
        Dab = self.sigma3 * Fab * self.kappaABij * ghs_eval
        dDab = Dab * dghs_rho/ghs_eval
        if type(a_hs[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')

        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])

        Dabij[self.indexabij] = Dab
        dDabij_drho[self.indexabij] = dDab
        
        KIJ = rho * (self.DIJ*Dabij)
        if Xass0 is None:
            sumD = np.sum(KIJ, axis = 1)
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)

        Xass = Xass_solver(self.nsites, KIJ, self.diagasso, Xass0)

        CIJ = rho * np.tile(Xass**2, (self.nsites, 1)).T * Dabij * self.DIJ
        CIJ[self.diagasso] += 1.
        dXass = dXass_drho(rho, Xass, self.DIJ, Dabij, dDabij_drho, CIJ)
        a[0] += np.dot(self.S, (np.log(Xass) - Xass/2 + 1/2))
        a[1] += np.dot(self.S, (1/Xass - 1/2) * dXass)
    else:
        Xass = None


    return a, Xass


def d2ares_drho(self, rho, temp_aux, Xass0=None):
    # Parameters needed for evaluating the helmothlz contributions
    beta, eps_beta, d, dia3, mes3, m2e2s3, deta_drho, Fab, betaxpmu2 = temp_aux
    eta = deta_drho * rho
    
    # Monomer contribution
    a_hs = d2ahs_deta(eta)

    # Chain contribution
    d2ghs = d2ghs_deta(eta)
    a_chain = d2achain_deta(self.ms, d2ghs)

    # Dispersion contribution
    a_disp = d2adisp_deta(self.ms, eta, dia3, self.ai, self.bi, m2e2s3, mes3)

    # Polar contribution
    apolar = np.zeros(3)
    if self.polar_bool:
        apolar = d2apolarJC_deta(eta, self.ms, dia3, betaxpmu2)

    # Total Helmolthz
    a = (self.ms*a_hs + a_chain) + a_disp + apolar

    a *= np.array([1, deta_drho, deta_drho**2]) 

    if self.assoc_bool:           
        ghs_eval = d2ghs[0]
        dghs_rho = d2ghs[1]*deta_drho
        d2ghs_rho = d2ghs[2]*deta_drho**2
        Dab = self.sigma3 * Fab * self.kappaABij * ghs_eval
        dDab = Dab * dghs_rho/ghs_eval
        d2Dab = Dab * d2ghs_rho/ghs_eval
        if type(a_hs[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            d2Dabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')

        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])
            d2Dabij_drho = np.zeros([self.nsites, self.nsites])

        Dabij[self.indexabij] = Dab
        dDabij_drho[self.indexabij] = dDab
        d2Dabij_drho[self.indexabij] = d2Dab

        KIJ = rho * (self.DIJ*Dabij)
        if Xass0 is None:
            sumD = np.sum(KIJ, axis = 1)
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)

        Xass = Xass_solver(self.nsites, KIJ, self.diagasso, Xass0)
        CIJ = rho * np.tile(Xass**2, (self.nsites, 1)).T * Dabij * self.DIJ
        CIJ[self.diagasso] += 1.
	
        dXass = dXass_drho(rho, Xass, self.DIJ, Dabij, dDabij_drho, CIJ)
        d2Xass = d2Xass_drho(rho, Xass, dXass, self.DIJ, Dabij, dDabij_drho,
                             d2Dabij_drho, CIJ)
        a[0] += np.dot(self.S, (np.log(Xass) - Xass/2 + 1/2))
        a[1] += np.dot(self.S, (1/Xass - 1/2) * dXass)
        a[2] += np.dot(self.S, - (dXass/Xass)**2 + d2Xass * (1/Xass - 1/2))
    else:
        Xass = None


    return a, Xass

def all_contributions(self, rho, temp_aux, Xass0=None):
    # Parameters needed for evaluating the helmothlz contributions
    beta, eps_beta, d, dia3, mes3, m2e2s3, deta_drho, Fab, betaxpmu2 = temp_aux
    eta = deta_drho * rho
    
    # Monomer contribution
    a_hs = d2ahs_deta(eta)
    aMono = self.ms*a_hs


    # Chain contribution
    d2ghs = d2ghs_deta(eta)
    aChain = d2achain_deta(self.ms, d2ghs)
       

    # Dispersion contribution
    aDisp = d2adisp_deta(self.ms, eta, dia3, self.ai, self.bi, m2e2s3, mes3)

    # Polar contribution
    aPolar = np.zeros(3)
    if self.polar_bool:
        aPolar = d2apolarJC_deta(eta, self.ms, dia3, betaxpmu2)

    # Total Helmolthz
    aux_vec = np.array([1, deta_drho, deta_drho**2]) 
    aMono *= aux_vec
    aChain *= aux_vec
    aDisp *= aux_vec
    aPolar *= aux_vec

    aAsso = aux_vec * 0.
    if self.assoc_bool:           
        ghs_eval = d2ghs[0]
        dghs_rho = d2ghs[1]*deta_drho
        d2ghs_rho = d2ghs[2]*deta_drho**2
        Dab = self.sigma3 * Fab * self.kappaABij * ghs_eval
        dDab = Dab * dghs_rho/ghs_eval
        d2Dab = Dab * d2ghs_rho/ghs_eval
        
        if type(a_hs[0]) == np.complex128:
            Dabij = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            dDabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')
            d2Dabij_drho = np.zeros([self.nsites, self.nsites], dtype = 'complex_')

        else:
            Dabij = np.zeros([self.nsites, self.nsites])
            dDabij_drho = np.zeros([self.nsites, self.nsites])
            d2Dabij_drho = np.zeros([self.nsites, self.nsites])
            
        Dabij[self.indexabij] = Dab
        dDabij_drho[self.indexabij] = dDab
        d2Dabij_drho[self.indexabij] = d2Dab

        KIJ = rho * (self.DIJ*Dabij)
        if Xass0 is None:
            sumD = np.sum(KIJ, axis = 1)
            Xass0 = (np.sqrt(1 + 4*sumD) - 1)/(2*sumD)

        Xass = Xass_solver(self.nsites, KIJ, self.diagasso, Xass0)

        CIJ = rho * np.tile(Xass**2, (self.nsites, 1)).T * Dabij * self.DIJ
        CIJ[self.diagasso] += 1.
        dXass = dXass_drho(rho, Xass, self.DIJ, Dabij, dDabij_drho, CIJ)
        d2Xass = d2Xass_drho(rho, Xass, dXass, self.DIJ, Dabij, dDabij_drho,
                             d2Dabij_drho, CIJ)
        aa1 = np.dot(self.S, (np.log(Xass) - Xass/2 + 1/2))
        aa2 = np.dot(self.S, (1/Xass - 1/2) * dXass)
        aa3 = np.dot(self.S, - (dXass/Xass)**2 + d2Xass * (1/Xass - 1/2))
        aAsso = np.array([aa1, aa2, aa3])

    else:
        Xass = None


    return aMono, aChain, aDisp,  aAsso, aPolar

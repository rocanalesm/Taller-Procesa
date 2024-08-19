from __future__ import division, print_function, absolute_import
import numpy as np

from .ares import ares, dares_drho, d2ares_drho, dares_dx, dares_dxrho, all_contributions

from .ideal.ideal import aideal, daideal_drho, d2aideal_drho, daideal_dx
from .ideal.ideal import daideal_dxrho

from .association.association_aux import association_config

from .routines.density_solver import density_newton
from ..constants import kb, Na, e, eps0
from ..pcsaft_pure import pcsaft_pure

R = Na * kb


class pcsaft_mix():
    '''
    PC-SAFT EoS for mixtures Object

    This object has implemented methods for phase equilibrium
    as for interfacial properties calculations.

    Parameters
    ----------
    mixture : object
        the mixture is created with mixture class

    Attributes
    ----------
    nc: integer
        number of components in the mixture
    ms: array_like
        number of chain segments
    sigma: array_like
        segment diameter  [m]
    eps: array_like
        dispersion energy [J]


    sigmaij: array_like
        segment diameter matrix [m]
    epsij: array_like
        dispersion energy matrix [J]

    eABij: array_like
        association energy matrix [J]
    kappaABij: array_like
        association volume
    sites: list
        triplet of number of association sites [B, P, N]

    mupol: array_like
        dipolar moment [Debye]
    xpol: array_like
        fraction of dipolar segment sites

    cii : array_like
        influence factor for SGT [J m^5 / mol^2]
    cij : array_like
        cross influence parameter matrix for SGT [J m^5 / mol^2]
    beta: array_like
        correction to cross-influence parameter matrix

    secondorder: bool
        bool to indicate if composition derivatives of fugacity coefficient
        are available
    secondordersgt: bool
        bool to indicate if composition derivatives of chemical potential
        are available


    Methods
    -------
    diameter : computes the diameter at a given temperature
    temperature_aux : computes temperature dependent parameters of the fluid
    ares: computes the residual dimensionless Helmholtz free energy
    dares_drho: computes the residual dimensionless Helmholtz free energy
                and its density first density derivative
    d2ares_drho: computes the residual dimensionless Helmholtz free energy
                 and its density first and second-density derivatives
    dares_dx: computes the residual dimensionless Helmholtz free energy and
              its composition derivatives
    dares_dxrho: computes the residual dimensionless Helmholtz free energy
                 and its composition and density derivatives
    afcn: computes total Helmholtz energy
    dafcn_drho: computes total Helmholtz energy and its density derivative
    d2afcn_drho: computes total Helmholtz energy and its density derivatives
    dafcn_dx: computes total Helmholtz energy and its composition derivative
    dafcn_dxrho: computes total Helmholtz energy and its composition and
                density derivative
    density: computes the density of the fluid
    pressure: computes the pressure
    dP_drho: computes pressure and its density derivative
    logfugmix: computes the fugacity coefficient of the mixture
    logfugef: computes the effective fugacity coefficients of the components
              in the mixture
    a0ad: computes adimensional Helmholtz density energy
    muad: computes adimensional chemical potential
    dmuad: computes the adimensional chemical potential and its derivatives
    dOm : computes adimensional Thermodynamic Grand Potential
    ci :  computes influence parameters matrix for SGT
    sgt_adim : computes adimensional factors for SGT
    beta_sgt: method for setting beta correction used in SGT
    EntropyR : computes the residual entropy of the fluid
    EnthalpyR : computes the residual enthalpy of the fluid
    CvR : computes the residual isochoric heat capacity
    CpR : computes the residual heat capacity
    speed_sound : computes the speed of sound

    Auxiliar methods (computed using temperature_aux output list)
    -------------------------------------------------------------
    density_aux : computes density
    afcn_aux : computes afcn
    dafcn_aux : computes dafcn_drho
    d2afcn_aux : computes d2afcn_drho
    pressure_aux : computes pressure
    dP_drho_aux : computes dP_drho
    logfugmix_aux : computes logfug
    a0ad_aux : compute a0ad
    muad_aux : computes muad
    dmuad_aux : computes dmuad
    dOm_aux : computes dOm
    
    Code template based on sgtpy python package (https://github.com/gustavochm/sgtpy)
    '''

    def __init__(self, mixture, compute_critical = True):

    
        self.mixture = mixture
        self.nc = mixture.nc
        self.Mw = np.asarray(mixture.Mw)
        self.diag_index = np.diag_indices(self.nc)
        
        # No-associative parameters
        self.ms = np.asarray(mixture.ms)
        self.sigmaT_bool = mixture.sigmaT_bool
        
        if not self.sigmaT_bool:
            self.sigma = np.asarray(mixture.sigma0)
        else:
            self.sigma = np.asarray(mixture.sigma0)
            self.sigma0 = np.asarray(mixture.sigma0)
            self.t1 = np.asarray(mixture.t1)
            self.t2 = np.asarray(mixture.t2)
            self.t3 = np.asarray(mixture.t3)
            self.t4 = np.asarray(mixture.t4)
        self.eps = np.asarray(mixture.eps)

        
        self.sigmaij0 = np.add.outer(self.sigma, self.sigma) / 2
        self.epsij0 = np.sqrt(np.outer(self.eps,self.eps))
        self.mij = np.outer(self.ms, self.ms)

        self.dxhi00_drho = np.pi / 6
        
        # Association parameters   
        self.eAB = np.asarray(mixture.eAB)
        self.kappaAB = np.asarray(mixture.kappaAB)
        self.eABij0 = np.add.outer(mixture.eAB, mixture.eAB)/2       
        sigma3 = self.sigma**3
        sigmaij3 = self.sigmaij0**3
        self.kappaAB = np.asarray(mixture.kappaAB)
        self.kappaABij = np.multiply.outer(sigma3, sigma3)
        self.kappaABij *= np.multiply.outer(self.kappaAB, self.kappaAB)
        self.kappaABij **= 0.5
        self.kappaABij /= sigmaij3   
        self.sitesmix = mixture.sitesmix
        S, DIJ, compindex, indexabij, indexab, nsites, \
        dxjdx, diagasso = association_config(self.sitesmix, self)
        assoc_bool = nsites != 0
        self.assoc_bool = assoc_bool
        if assoc_bool:
            self.S = S
            self.DIJ = DIJ
            self.compindex = compindex
            self.indexab = indexab
            self.indexabij = indexabij
            self.nsites = nsites
            self.dxjdx = dxjdx
            self.diagasso = diagasso

            
        # Polar parameters
        self.mupol = np.asarray(mixture.mupol)
        self.xpol = np.asarray(mixture.xpol)
        self.polar_bool = np.any(self.xpol != 0)
        if self.polar_bool:
            # 1 D = 3.33564e-30 C * m
            # 1 C^2 = 9e9 N m^2
            cte = (3.33564e-30)**2 * (9e9)
            self.mxpmupol2 = self.ms * self.xpol * self.mupol**2 * cte
        else:
            self.mxpmupol2 = 0
            
        # Electrolyte parameters 
        self.z = np.asarray(mixture.z)
        self.elec_bool = np.any(self.z != 0)
        self.index_ions = (self.z != 0)
        self.index_solv = (self.z == 0)
        self.Mw_solv = self.Mw[self.index_solv]
        cation_cation = np.outer(self.z > 0, self.z > 0)
        anion_anion = np.outer(self.z < 0, self.z < 0)
        self.index_disp = (cation_cation + anion_anion) == 0
        self.index_disp += np.outer(self.ms != 1, self.ms != 1)

        
        # eri = er0 + er1*T
        self.er0 = np.asarray(mixture.er0)
        self.er1 = np.asarray(mixture.er1)
        self.er2 = np.asarray(mixture.er2)
        self.e2_eps0 = e**2  / eps0
            
        # Binary parameters
        # kij matrix: kij = kij0 + kij1 T + kij2 T^2 + kij3 / T
        self.KIJ0saft = mixture.KIJ0saft
        self.KIJ1saft = mixture.KIJ1saft
        self.KIJ2saft = mixture.KIJ2saft
        self.KIJ3saft = mixture.KIJ3saft
        
        # kepsij matrix: kepsij = kepsij0 + kepsij1 T + kepsij2 T^2 + kepsij3 / T
        self.KepsIJ0saft = mixture.KepsIJ0saft
        self.KepsIJ1saft = mixture.KepsIJ1saft
        self.KepsIJ2saft = mixture.KepsIJ2saft
        self.KepsIJ3saft = mixture.KepsIJ3saft
        
        # lij matrix: lij = lij0 + lij1 T + lij2 T^2 + lij3 / T
        self.LIJ0saft = mixture.LIJ0saft
        self.LIJ1saft = mixture.LIJ1saft
        self.LIJ2saft = mixture.LIJ2saft
        self.LIJ3saft = mixture.LIJ3saft


        # DGT parameters
        self.cii = mixture.cii
        self.secondorder = False
        self.secondordersgt = True
        self.beta0 = np.zeros([self.nc, self.nc])
        self.beta1 = np.zeros([self.nc, self.nc])
        self.beta2 = np.zeros([self.nc, self.nc])
        self.beta3 = np.zeros([self.nc, self.nc])

        # entropy scaling calculation
        self.viscosity_parameters = np.array(mixture.viscosity_parameters)
        
        # creating pure fluid's eos
        pure_eos = []
        for i, component in enumerate(self.mixture.components):
            model = pcsaft_pure(component, compute_critical=compute_critical)
            pure_eos.append(model)
        self.pure_eos = pure_eos


    def temperature_aux(self, T):
        """
        temperature_aux(T)

        Method that computes temperature-dependent parameters.

        Parameters
        ----------

        T : float
            Absolute temperature [K]

        Returns
        -------
        temp_aux : list
             list of computed parameters
        """
        if self.sigmaT_bool:
            self.sigma = T**0 * self.sigma0
            self.sigma += self.t1*np.exp(self.t2*T)
            self.sigma += self.t3*np.exp(self.t4*T)
            self.sigmaij0 = np.add.outer(self.sigma, self.sigma) / 2
            sigma3 = self.sigma**3
            sigmaij3 = self.sigmaij0**3
            self.kappaABij = np.multiply.outer(sigma3, sigma3)
            self.kappaABij *= np.multiply.outer(self.kappaAB, self.kappaAB)
            self.kappaABij **= 0.5
            self.kappaABij /= sigmaij3
            
        # computing temperature dependent kij and epsij
        kij = self.KIJ0saft + self.KIJ1saft*T + self.KIJ2saft*T**2
        kij += self.KIJ3saft/T
        self.epsij = self.epsij0 * (1. - kij)
        
        # computing temperature dependent kepsij and eABij
        kepsij = self.KepsIJ0saft + self.KepsIJ1saft*T + self.KepsIJ2saft*T**2
        kepsij += self.KepsIJ3saft/T
        self.eABij = self.eABij0 * (1. - kepsij)    
        
        # computing temperature dependent lij and sigmaij
        lij = self.LIJ0saft + self.LIJ1saft*T + self.LIJ2saft*T**2
        lij += self.LIJ3saft/T
        self.sigmaij = self.sigmaij0 * (1. - lij)
        sigmaij3 = self.sigmaij**3
        

        beta = 1 / (kb * T)
        eps_beta = self.eps*beta

        # Hard sphere auxiliary parameters
        dii = self.sigma*(1 - 0.12*np.exp(-3*eps_beta*np.logical_or(self.z==0, self.ms != 1)))
        di03 = np.power.outer(dii, np.arange(4))
        aux_dii = np.multiply.outer(dii, dii)/np.add.outer(dii, dii)
        aux_dii2 = aux_dii**2
        
        # Dispersion auxiliary parameters
        epsij_beta = self.epsij*beta
        es3ij = epsij_beta * sigmaij3 * self.index_disp
        e2s3ij = epsij_beta**2 * sigmaij3 * self.index_disp
        mes3ij = self.mij * es3ij 
        me2s3ij = self.mij * e2s3ij 
        
        # Association auxiliary parameters
        Fab = np.exp(beta * self.eABij) - 1.
        aux_Dab = sigmaij3 * Fab * self.kappaABij
        
        # Polar auxiliary parameters
        aux1_polarij = np.zeros_like(aux_Dab)
        aux2_polarijk = np.zeros([self.nc, self.nc, self.nc])
        dij3 = np.zeros_like(aux_Dab)
        dijk3 = np.zeros_like(aux2_polarijk)
        if self.polar_bool:
            mxpmu2 = beta * self.mxpmupol2
            dij = np.add.outer(dii, dii) / 2
            #dij *= (1. - lij)
            dij3 = dij**3
            
            mult = np.multiply.outer(dij, dij)
            listpolar = np.arange(self.nc)
            dijk3 = mult[listpolar, :, listpolar] * dij

            aux1_polarij = np.multiply.outer(mxpmu2, mxpmu2)
            aux2_polarijk = np.multiply.outer(mxpmu2, aux1_polarij)
            
        # Electrolite auxiliary parameters
        betae2_e0 = 0.
        if self.elec_bool:
            self.eri = self.er0 + self.er1 * T + self.er2 * np.log(T)
            betae2_e0 = self.e2_eps0 * beta
            self.eri_ions = self.eri[self.index_ions]
            self.eri_solv = self.eri[self.index_solv]

        
        temp_aux = [beta, di03, aux_dii, aux_dii2, es3ij, e2s3ij, 
                    mes3ij, me2s3ij, aux_Dab, aux1_polarij, aux2_polarijk, 
                    dij3, dijk3, dii, betae2_e0]
        return temp_aux

    def ares(self, x, rho, T, Xass0=None):
        """
        ares(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the mixture.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: float
           residual dimensionless Helmholtz free energy [Adim]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = ares(self, x, rho, temp_aux, Xass0)
        return a

    def dares_drho(self, x, rho, T, Xass0=None):
        """
        dares_drho(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the mixture
        and its first density derivative.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           residual dimensionless Helmholtz free energy [Adim, m^3]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = dares_drho(self, x, rho, temp_aux, Xass0)
        return a

    def d2ares_drho(self, x, rho, T, Xass0=None):
        """
        d2ares_drho(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the mixture
        and its first and second-density derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           residual dimensionless Helmholtz free energy [Adim, m^3, m^6]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = d2ares_drho(self, x, rho, temp_aux, Xass0)
        return a

    def dares_dx(self, x, rho, T, Xass0=None):
        """
        dares_dx(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the mixture
        and its composition derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: float
           residual dimensionless Helmholtz free energy [Adim]
        ax: array_like
           composition derivatives of residual dimensionless Helmholtz
           free energy [Adim]
        """
        temp_aux = self.temperature_aux(T)
        a, ax, Xass = dares_dx(self, x, rho, temp_aux, Xass0)
        return a, ax

    def dares_dxrho(self, x, rho, T, Xass0=None):
        """
        dares_dx(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the mixture
        and its density and composition derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           residual dimensionless Helmholtz free energy [Adim, m^3]
        ax: array_like
           composition derivatives of residual dimensionless Helmholtz
           free energy [Adim]
        """
        temp_aux = self.temperature_aux(T)
        a, ax, Xass = dares_dxrho(self, x, rho, temp_aux, Xass0)
        return a, ax

    def afcn_aux(self, x, rho, temp_aux, Xass0=None):
        """
        afcn_aux(x, rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the mixture.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature-dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: float
           total Helmholtz free energy [J/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        a, Xass = ares(self, x, rho, temp_aux, Xass0)
        a += aideal(x, rho, beta)
        a *= (Na/beta)
        return a, Xass

    def dafcn_drho_aux(self, x, rho, temp_aux, Xass0=None):
        """
        dafcn_drho_aux(x, rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its first density derivative.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature-dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           toal Helmholtz free energy [J/mol, J m^3/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        a, Xass = dares_drho(self, x, rho, temp_aux, Xass0)
        a += daideal_drho(x, rho, beta)
        a *= (Na/beta)
        return a, Xass

    def d2afcn_drho_aux(self, x, rho, temp_aux, Xass0=None):
        """
        d2afcn_drho_aux(x, rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its first and second-density derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           toal Helmholtz free energy [J/mol, J m^3/mol, J m^6/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        a, Xass = d2ares_drho(self, x, rho, temp_aux, Xass0)
        a += d2aideal_drho(x, rho, beta)
        a *= (Na/beta)
        return a, Xass

    def dafcn_dx_aux(self, x, rho, temp_aux, Xass0=None):
        """
        dafcn_dx_aux(x, rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its composition derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: float
           total dimensionless Helmholtz free energy [J/mol]
        ax: array_like
           composition derivatives of total dimensionless Helmholtz
           free energy [J/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        ar, aresx, Xass = dares_dx(self, x, rho, temp_aux, Xass0)
        aideal, aidealx = daideal_dx(x, rho, beta)
        a = (ar + aideal)
        a *= (Na/beta)
        ax = (aresx + aidealx)
        ax *= (Na/beta)
        return a, ax, Xass

    def dafcn_dxrho_aux(self, x, rho, temp_aux, Xass0=None):
        """
        dafcn_dxrho_aux(x, rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its composition and density derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           total dimensionless Helmholtz free energy [J/mol, J m^3/mol]
        ax: array_like
           composition derivatives of total dimensionless Helmholtz
           free energy [J/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        ar, aresx, Xass = dares_dxrho(self, x, rho, temp_aux, Xass0)
        aideal, aidealx = daideal_dxrho(x, rho, beta)
        a = (ar + aideal)
        a *= (Na/beta)
        ax = (aresx + aidealx)
        ax *= (Na/beta)
        return a, ax, Xass

    def afcn(self, x, rho, T, Xass0=None):
        """
        afcn(x, rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the mixture.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: float
           total Helmholtz free energy [J/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = self.afcn_aux(x, rho, temp_aux, Xass0)
        return a

    def dafcn_drho(self, x, rho, T, Xass0=None):
        """
        dafcn_drho(x, rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its first density derivative.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           toal Helmholtz free energy [J/mol, J m^3/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = self.dafcn_drho_aux(x, rho, temp_aux, Xass0)
        return a

    def d2afcn_drho(self, x, rho, T, Xass0=None):
        """
        d2afcn_drho(x, rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its first and second-density derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           toal Helmholtz free energy [J/mol, J m^3/mol, J m^6/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = self.d2afcn_drho_aux(x, rho, temp_aux, Xass0)
        return a

    def dafcn_dx(self, x, rho, T, Xass0=None):
        """
        dafcn_dx(x, rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its composition derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: float
           total dimensionless Helmholtz free energy [J/mol]
        ax: array_like
           composition derivatives of total dimensionless Helmholtz
           free energy [J/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, ax, Xass = self.dafcn_dx_aux(x, rho, temp_aux, Xass0)
        return a, ax

    def dafcn_dxrho(self, x, rho, T, Xass0=None):
        """
        dafcn_dxrho(x, rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the mixture
        and its composition and density derivatives.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a: array_like
           total dimensionless Helmholtz free energy [J/mol, J m^3/mol]
        ax: array_like
           composition derivatives of total dimensionless Helmholtz
           free energy [J/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, ax, Xass = self.dafcn_dxrho_aux(x, rho, temp_aux, Xass0)
        return a, ax

    def density_aux(self, x, temp_aux, P, state, rho0=None, Xass0=None):
        """
        density_aux(x, temp_aux, P, state)
        Method that computes the density of the mixture at given composition,
        temperature, pressure and aggregation state.

        Parameters
        ----------
        x: array_like
            molar fraction array
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapor phase
        rho0 : float, optional
            initial guess to compute density root [mol/m^3]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        density: float
            density [mol/m^3]
        Xass : array
            computed fraction of nonbonded sites
        """
        etamax = 0.7404804896930609
        rhomax = (6 * etamax) 
        rhomax /= np.dot(x, (self.ms * np.pi * self.sigma**3))
        rhomax /= Na
                
        if rho0 is None:
            if state == 'L':
                eta0 = 0.5 * P**0
                rho0 = rhomax * eta0 / etamax
            elif state == 'V':
                beta = temp_aux[0]
                rho0 = P*beta/Na

        rho, Xass = density_newton(rho0, x, temp_aux, P, Xass0, rhomax.real, state, self)

        return rho, Xass

    def density(self, x, T, P, state, rho0=None, Xass0=None):
        """
        density(x, T, P, state)
        Method that computes the density of the mixture at a given composition,
        temperature, pressure, and aggregation state.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapor phase
        rho0 : float, optional
            initial guess to compute density root [mol/m^3]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        density: float
            density [mol/m^3]
        """
        temp_aux = self.temperature_aux(T)
        rho, Xass = self.density_aux(x, temp_aux, P, state, rho0, Xass0)
        return rho

    def pressure_aux(self, x, rho, temp_aux, Xass0=None):
        """
        pressure_aux(x, rho, temp_aux, Xass0)

        A method that computes the pressure at a given composition,
        density [mol/m3] and temperature [K]

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            density [mol/m3]
        temp_aux : list
            temperature-dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        Xass : array
            computed fraction of nonbonded sites
        """
        rhomolecular = Na * rho
        da, Xass = self.dafcn_drho_aux(x, rhomolecular, temp_aux, Xass0)
        afcn, dafcn = da
        Psaft = rhomolecular**2 * dafcn / Na
        return Psaft, Xass

    def dP_drho_aux(self, x, rho, temp_aux, Xass0=None):
        """
        dP_drho_aux(rho, temp_aux, Xass0)

        Method that computes the pressure and its density derivative at a given
        composition, density [mol/m3], and temperature [K]

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            density [mol/m3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        dP: float
            derivate of pressure respect density [Pa m^3 / mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        rhomolecular = Na * rho
        da, Xass = self.d2afcn_drho_aux(x, rhomolecular, temp_aux, Xass0)
        afcn, dafcn, d2afcn = da
        Psaft = rhomolecular**2 * dafcn / Na
        dPsaft = 2 * rhomolecular * dafcn + rhomolecular**2 * d2afcn
        return Psaft, dPsaft, Xass

    def pressure(self, x, rho, T, Xass0=None):
        """
        pressure(x, rho, T, Xass0)

        Method that computes the pressure at a given composition,
        density [mol/m3] and temperature [K]

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            density [mol/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        """
        temp_aux = self.temperature_aux(T)
        Psaft, Xass = self.pressure_aux(x, rho, temp_aux, Xass0)
        return Psaft

    def dP_drho(self, x, rho, T, Xass0=None):
        """
        dP_drho(rho, T, Xass0)

        Method that computes the pressure and its density derivative at a given
        composition, density [mol/m3], and temperature [K]

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            density [mol/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        dP: float
            derivate of pressure respect density [Pa m^3 / mol]
        """
        temp_aux = self.temperature_aux(T)
        Psaft, dPsaft, Xass = self.dP_drho_aux(x, rho, temp_aux, Xass0)
        return Psaft, dPsaft

    def logfugmix_aux(self, x, temp_aux, P, state, v0=None, Xass0=None):
        """
        logfugmix_aux(x, temp_aux, P, state, v0, Xass0)

        Method that computes the fugacity coefficient of the mixture at a given
        composition, temperature, and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        temp_aux : list
            temperature-dependent parameters computed with temperature_aux(T)
        P: float
            pressure [Pa]
        state: string
            'L' for the liquid phase and 'V' for vapor phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        lnphi: float
            fugacity coefficient of the mixture
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

        v = 1./rho
        rhomolecular = Na * rho
        ar, Xass = ares(self, x, rhomolecular, temp_aux, Xass)
        Z = P * v / RT
        lnphi = ar + (Z - 1.) - np.log(Z)
        return lnphi, v, Xass

    def logfugmix(self, x, T, P, state, v0=None, Xass0=None):
        """
        logfugmix(x, T, P, state, v0, Xass0)

        Method that computes the fugacity coefficient of the mixture at a given
        composition, temperature, and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T: float
            absolute temperature [K]
        P: float
            pressure [Pa]
        state: string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        lnphi: float
            fugacity coefficient of the mixture
        v: float
            the computed volume of the phase [m^3/mol]
        """
        temp_aux = self.temperature_aux(T)
        lnphi, v, Xass = self.logfugmix_aux(x, temp_aux, P, state, v0, Xass0)
        return lnphi, v

    def logfugef_aux(self, x, temp_aux, P, state, v0=None, Xass0=None):
        """
        logfugef_aux(x, temp_aux, P, state, v0, Xass0)

        A method that computes the effective fugacity coefficient of the
        components in the mixture at a given composition, temperature
        and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        temp_aux : list
            temperature-dependent parameters computed with temperature_aux(T)
        P: float
            pressure [Pa]
        state: string
            'L' for the liquid phase and 'V' for the vapor phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        lnphi: array_like
            effective fugacity coefficient of the components
        v: float
            the computed volume of the phase [m^3/mol]
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
        v = 1./rho
        rhomolecular = Na * rho
        ar, daresx, Xass = dares_dx(self, x, rhomolecular, temp_aux, Xass)
        Z = P * v / RT
        mures = ar + (Z - 1.) + daresx - np.dot(x, daresx)
        lnphi = mures - np.log(Z)
        return lnphi, v, Xass

    def logfugef(self, x, T, P, state, v0=None, Xass0=None):
        """
        logfugef(x, T, P, state, v0, Xass0)

        A method that computes the effective fugacity coefficient of the
        components in the mixture at a given composition, temperature
        and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T: float
            absolute temperature [K]
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
        lnphi: array_like
            effective fugacity coefficient of the components
        v: float
            the computed volume of the phase [m^3/mol]
        """
        temp_aux = self.temperature_aux(T)
        lnphi, v, Xass = self.logfugef_aux(x, temp_aux, P, state, v0, Xass0)
        return lnphi, v

    def a0ad_aux(self, rhoi, temp_aux, Xass0=None):
        """
        a0ad_aux(rhoi, temp_aux, Xass0)

        A method that computes the Helmholtz density energy divided by RT at
        given density vector and temperature.

        Parameters
        ----------

        rhoi : array_like
            density vector [mol/m^3]
        temp_aux : list
            temperature-dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        a0ad: float
            Helmholtz density energy divided by RT [mol/m^3]
        Xass : array
            computed fraction of nonbonded sites
        """

        rho = np.sum(rhoi)
        x = rhoi/rho
        rhomolecular = Na * rho
        beta = temp_aux[0]

        a, Xass = ares(self, x, rhomolecular, temp_aux, Xass0)
        a += aideal(x, rhomolecular, beta)

        a0 = a*rho

        return a0, Xass

    def a0ad(self, rhoi, T, Xass0=None):
        """
        a0ad(rhoi, T, Xass0)

        Method that computes the Helmholtz density energy divided by RT at
        given density vector and temperature.

        Parameters
        ----------

        rhoi : array_like
            density vector [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a0ad: float
            Helmholtz density energy divided by RT [mol/m^3]
        """
        temp_aux = self.temperature_aux(T)
        a0, Xass = self.a0ad_aux(rhoi, temp_aux, Xass0)
        return a0

    def muad_aux(self, rhoi, temp_aux, Xass0=None):
        """
        muad_aux(rhoi, temp_aux, Xass0)

        A method that computes the dimensionless chemical potential at a given
        density vector and temperature.

        Parameters
        ----------
        rhoi : array_like
            density vector [mol/m^3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        muad: array_like
            chemical potential [Adim]
        Xass : array
            computed fraction of nonbonded sites
        """
        rho = np.sum(rhoi)
        x = rhoi/rho
        rhom = Na * rho

        beta = temp_aux[0]
        ares, aresx, Xass = dares_dxrho(self, x, rhom, temp_aux, Xass0)
        aideal, aidealx = daideal_dxrho(x, rhom, beta)
        afcn, dafcn = (ares + aideal)
        ax = (aresx + aidealx)
        Z = dafcn * rhom
        mu = afcn + ax - np.dot(x, ax) + (Z)
        return mu, Xass

    def muad(self, rhoi, T, Xass0=None):
        """
        muad(rhoi, T, Xass0)

        A method that computes the dimensionless chemical potential at a given
        density vector and temperature.

        Parameters
        ----------
        rhoi : array_like
            density vector [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        muad: array_like
            chemical potential [Adim]
        """
        temp_aux = self.temperature_aux(T)
        mu, Xass = self.muad_aux(rhoi, temp_aux, Xass0)
        return mu

    def dmuad_aux(self, rhoi, temp_aux, Xass0=None):
        """
        dmuad_aux(rhoi, temp_aux, Xass0)

        A method that computes the chemical potential and its numerical
        derivative at a given density vector and temperature.

        Parameters
        ----------
        rhoi : array_like
            density vector [mol/m^3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        muad: array_like
            chemical potential [J/mol]
        dmuad: array_like
            derivatives of the chemical potential with respect to rhoi [J m^3/mol^2]
        Xass : array
            computed fraction of nonbonded sites
        """
        nc = self.nc
        h = np.finfo(1.0).eps
        diff = h * np.eye(nc) * 1j 

        mu, Xass = self.muad_aux(rhoi, temp_aux, Xass0)

        if isinstance(Xass, np.ndarray):
            Xass = Xass.astype('complex128')
            
        arr = []
        for i in range(nc):
            muad1, _ = self.muad_aux(rhoi.astype('complex128') + diff[i], temp_aux, Xass)
            arr.append(muad1.imag/h)

        dmu = np.column_stack(arr)
        return mu, dmu, Xass

    def dmuad(self, rhoi, T, Xass0=None):
        """
        dmuad(rhoi, T, Xass0)

        A method that computes the chemical potential and its numerical
        derivative at a given density vector and temperature.

        Parameters
        ----------
        rhoi : array_like
            density vector [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        muad: array_like
            chemical potential [J/mol]
        dmuad: array_like
            derivatives of the chemical potential with respect to rhoi [J m^3/mol^2]
        """
        temp_aux = self.temperature_aux(T)
        mu, dmu, Xass = self.dmuad_aux(rhoi, temp_aux, Xass0)
        return mu, dmu

    def dOm_aux(self, rhoi, temp_aux, mu, Psat, Xass0=None):
        """
        dOm_aux(rhoi, temp_aux, mu, Psat, Xass0)

        A method that computes the Thermodynamic Grand potential
        at given density and temperature.

        Parameters
        ----------
        rhoi : array_like
            density vector [mol/m^3]
        temp_aux : list
            temperature dependent parameters computed with temperature_aux(T)
        mu : float
            adimensional chemical potential at equilibrium
        Psat : float
            equilibrium pressure divided by RT [Pa mol / J]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        dom: float
            Thermodynamic Grand potential [Pa mol / J]
        Xass : array
            computed fraction of nonbonded sites
        """
        a0ad, Xass = self.a0ad_aux(rhoi, temp_aux, Xass0)
        dom = a0ad - np.sum(np.nan_to_num(rhoi*mu)) + Psat
        return dom, Xass

    def dOm(self, rhoi, T, mu, Psat, Xass0=None):
        """
        dOm(rhoi, T, mu, Psat, Xass0)

        A method that computes the Thermodynamic Grand potential
        at given density and temperature.

        Parameters
        ----------
        rhoi : array_like
            density vector [mol/m^3]
        T : float
            absolute temperature [K]
        mu : float
            adimensional chemical potential at equilibrium
        Psat : float
            equilibrium pressure divided by RT [Pa mol / J]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        Out: float
            Thermodynamic Grand Potential [Pa]
        """
        temp_aux = self.temperature_aux(T)
        dom, Xass = self.dOm_aux(rhoi, temp_aux, mu, Psat, Xass0)
        return dom

    def sgt_adim(self, T):
        '''
        sgt_adim(T)

        A method that evaluates adimensional factor for temperature, pressure,
        density, tension, and distance for interfacial properties computations
        with SGT.

        Parameters
        ----------
        T : float
        absolute temperature [K]

        Returns
        -------
        Tfactor : float
            factor to obtain dimensionless temperature (K -> K)
        Pfactor : float
            factor to obtain dimensionless pressure (Pa -> Pa/RT)
        rofactor : float
            factor to obtain dimensionless density (mol/m3 -> mol/m3)
        tenfactor : float
            factor to obtain dimensionless surface tension (mN/m)
        zfactor : float
            factor to obtain dimensionless distance  (Amstrong -> m)
        '''
        beta = 1 / (kb*T)
        RT = (Na/beta)
        cii0 = np.polyval(self.cii[0], T)  # computing first component cii

        Tfactor = 1.
        Pfactor = 1. / RT
        rofactor = 1.
        tenfactor = np.sqrt(cii0*RT) * 1000  # To give tension in mN/m
        zfactor = 10**-10 * np.sqrt(RT / cii0)
        return Tfactor, Pfactor, rofactor, tenfactor, zfactor

    def beta_sgt(self, beta0, beta1=None, beta2=None, beta3=None):
        r"""
        beta_sgt(beta)

        A method that adds a beta correction for cross-influence parameters is used
        in SGT. The beta correction is computed as follows:

        .. math::
            \beta_{ij} =  \beta_{ij,0} + \beta_{ij,1} \cdot T +  \beta_{ij,2} \cdot T^2 + \frac{\beta_{ij,3}}{T}

        Parameters
        ----------
        beta0 : array_like
            beta0 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [Adim]
        beta1 : array_like, optional
            beta1 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K].
            If None, then a zero matrix is assumed.
        beta2 : array_like, optional
            beta2 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K^2].
            If None, then a zero matrix is assumed.
        beta3 : array_like, optional
            beta3 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [K]
            If None, then a zero matrix is assumed.

        """
        nc = self.nc

        Beta0 = np.asarray(beta0)
        shape = Beta0.shape
        isSquare = shape == (nc, nc)
        isSymmetric = np.allclose(Beta0, Beta0.T)
        diagZero = np.all(np.diagonal(Beta0) == 0.)
        if isSquare and isSymmetric and diagZero:
            self.beta0 = Beta0
        else:
            raise Exception('beta0 matrix is not square, symmetric or diagonal==0')

        if beta1 is None:
            Beta1 = np.zeros([nc, nc])
            self.beta1 = Beta1
        else:
            Beta1 = np.asarray(beta1)
            shape = Beta1.shape
            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(Beta1, Beta1.T)
            diagZero = np.all(np.diagonal(Beta1) == 0.)
            if isSquare and isSymmetric and diagZero:
                self.beta1 = Beta1
            else:
                raise Exception('beta1 matrix is not square, symmetric or diagonal==0')
        if beta2 is None:
            Beta2 = np.zeros([nc, nc])
            self.beta2 = Beta2
        else:
            Beta2 = np.asarray(beta2)
            shape = Beta2.shape
            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(Beta2, Beta2.T)
            diagZero = np.all(np.diagonal(Beta2) == 0.)
            if isSquare and isSymmetric and diagZero:
                self.beta2 = Beta2
            else:
                raise Exception('beta2 matrix is not square, symmetric or diagonal==0')

        if beta3 is None:
            Beta3 = np.zeros([nc, nc])
            self.beta3 = Beta3
        else:
            Beta3 = np.asarray(beta3)
            shape = Beta3.shape
            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(Beta3, Beta3.T)
            diagZero = np.all(np.diagonal(Beta3) == 0.)
            if isSquare and isSymmetric and diagZero:
                self.beta3 = Beta3
            else:
                raise Exception('beta3 matrix is not square, symmetric or diagonal==0')

    def set_betaijsgt(self, i, j, beta0, beta1=0., beta2=0., beta3=0.):
        r"""
        set_betaijsgt(i,j, beta0, beta1, beta2, beta3)

        A method that sets the betaij correction for cross-influence parameter
        between components i and j.
        The beta correction is computed as follows:

        .. math::
            \beta_{ij} =  \beta_{ij,0} + \beta_{ij,1} \cdot T +  \beta_{ij,2} \cdot T^2 + \frac{\beta_{ij,3}}{T}

        Parameters
        ----------
        i : int
            index of component i.
        j : int
            index of component j.
        beta0 : float
            beta0 value between component i and j [Adim]
        beta1 : float, optional
            beta1 value between component i and j [1/K]. Default to zero.
        beta2 : float, optional
            beta2 value between component i and j [1/K^2]. Default to zero.
        beta3 : float, optional
            beta3 value between component i and j [K]. Default to zero.

        """
        typei = type(i) == int
        typej = type(j) == int

        nc = self.nc
        nc_i = 0 <= i <= (nc - 1)
        nc_j = 0 <= j <= (nc - 1)

        i_j = i != j

        if (not nc_i) or (not nc_j):
            raise Exception('Index i or j bigger than (nc-1)')
        if not i_j:
            raise Exception('Cannot set betaij for i=j')

        if typei and typej and nc_i and nc_j and i_j:
            self.beta0[i, j] = beta0
            self.beta0[j, i] = beta0

            self.beta1[i, j] = beta1
            self.beta1[j, i] = beta1

            self.beta2[i, j] = beta2
            self.beta2[j, i] = beta2

            self.beta3[i, j] = beta3
            self.beta3[j, i] = beta3

    def ci(self, T):
        """
        A method that computes the matrix of cij interaction parameter for SGT at
        given temperature.

        Parameters
        ----------
        T : float
            absolute temperature [K]

        Returns
        -------
        ci : array_like
            influence parameter matrix at a given temperature [J m^5 / mol^2]

        """
        n = self.nc
        ci = np.zeros(n)
        for i in range(n):
            ci[i] = np.polyval(self.cii[i], T)
        cij = np.sqrt(np.outer(ci, ci))

        beta = self.beta0 + self.beta1*T + self.beta2*T**2 + self.beta3/T

        cij *= (1 - beta)
        return cij

    def EntropyR(self, x, T, P, state, v0=None, Xass0=None):
        """
        EntropyR(x, T, P, state, v0, Xass0, T_step)

        Method that computes the residual entropy (NPT) of the mixture at a given
        temperature and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for the liquid phase and 'V' for the vapor phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites


        Returns
        -------
        Sr : float
            residual entropy [J/mol K]

        """

        temp_aux = self.temperature_aux(T)
        if v0 is None:
            rho0 = None
        else:
            rho0 = 1./v0

        rho, Xass = self.density_aux(x, temp_aux, P, state, rho0, Xass0)

        v = 1/rho
        rhomolecular = Na * rho
        a, Xass = ares(self, x, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta
        Z = P * v / RT

        if type(Xass) == np.ndarray:
            Xass = Xass.astype('complex128')
        h = np.finfo(1.0).eps
        temp_aux_h = self.temperature_aux(T + h * 1j )

        ah, _ = ares(self, x.astype('complex128'), rhomolecular.astype('complex128'), temp_aux_h, Xass)

        dadT = ah.imag/h
        Sr_TVN = -T*dadT - a  # residual entropy (TVN) divided by R
        Sr_TPN = Sr_TVN + np.log(Z)  # residual entropy (TPN) divided by R
        Sr_TPN *= R
        return Sr_TPN

    def EntropyR_NrhoT(self, x, rho, T, Xass0 = None):
        """
        EntropyR_NrhoT(x, x, rho, T, Xass0 = None)

        A method that computes the residual entropy (NVT) at a given density and
        temperature.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho: float
            density [mol/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites


        Returns
        -------
        Sr : float
            residual entropy [J/mol K]

        """

        temp_aux = self.temperature_aux(T)
        rhomolecular = Na * rho
        a, Xass = ares(self, x, rhomolecular, temp_aux, Xass0)

        if isinstance(Xass, np.ndarray):
            Xass = Xass.astype('complex128')
        h = np.finfo(1.0).eps
        temp_aux_h = self.temperature_aux(T + h * 1j )

        ah, _ = ares(self, x.astype('complex128'), rhomolecular.astype('complex128'), temp_aux_h, Xass)

        dadT = ah.imag/h
        Sr_TVN = -T*dadT - a  # residual entropy (TVN) divided by R
        Sr_TVN *= R
        return Sr_TVN

    def EnthalpyR(self, x, T, P, state, v0=None, Xass0=None):
        """
        EnthalpyR(x, T, P, state, v0, Xass0, T_step)

        Method that computes the residual enthalpy of the mixture at a given
        temperature and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for the liquid phase and 'V' for the vapor phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites

        Returns
        -------
        Hr : float
            residual enthalpy [J/mol]

        """
        temp_aux = self.temperature_aux(T)
        if v0 is None:
            rho0 = None
        else:
            rho0 = 1./v0
        rho, Xass = self.density_aux(x, temp_aux, P, state, rho0, Xass0)

        v = 1/rho
        rhomolecular = Na * rho
        a, Xass = ares(self, x, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta
        Z = P * v / RT

        if isinstance(Xass, np.ndarray):
            Xass = Xass.astype('complex128')
        h = np.finfo(1.0).eps
        temp_aux_h = self.temperature_aux(T + h * 1j )

        ah, _ = ares(self, x, rhomolecular.astype('complex128'), temp_aux_h, Xass)

        dadT = ah.imag/h

        Sr_TVN = -T*dadT - a          # residual entropy (TVN) divided by R
        Hr_TPN = a + Sr_TVN + Z - 1.  # residual enthalpy divided by RT
        Hr_TPN *= RT
        return Hr_TPN

    def CvR(self, x, rho, T, Xass0=None, T_step=0.1):
        """
        CvR(x, rho, T, Xass0, T_step)

        Method that computes the residual isochoric heat capacity of the
        mixture at a given density and temperature.

        Parameters
        ----------
        x: array_like
            molar fraction array
        rho : float
            density [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites
        T_step: float, optional
            Step to compute temperature numerical derivates of Helmholtz
            free energy

        Returns
        -------
        Cv: float
            isochoric heat capacity [J/mol K]
        """
        temp_aux = self.temperature_aux(T)

        rhomolecular = Na * rho

        a, Xass = ares(self, x, rhomolecular, temp_aux, Xass0)

        h = T_step
        temp_aux1 = self.temperature_aux(T+h)
        temp_aux2 = self.temperature_aux(T+2*h)
        temp_aux_1 = self.temperature_aux(T-h)
        temp_aux_2 = self.temperature_aux(T-2*h)

        a1, Xass1 = ares(self, x, rhomolecular, temp_aux1, Xass)
        a2, Xass2 = ares(self, x, rhomolecular, temp_aux2, Xass)
        a_1, Xass_1 = ares(self, x, rhomolecular, temp_aux_1, Xass)
        a_2, Xass_2 = ares(self, x, rhomolecular, temp_aux_2, Xass)

        dFdT = (a_2/12 - 2*a_1/3 + 2*a1/3 - a2/12)/h
        d2FdT = (-a_2/12 + 4*a_1/3 - 5*a/2 + 4*a1/3 - a2/12)/h**2

        Cvr_TVN = -T**2*d2FdT - 2*T*dFdT  # residual isochoric heat capacity
        Cvr_TVN *= R
        return Cvr_TVN

    def CpR(self, x, T, P, state, v0=None, Xass0=None, T_step=0.1):
        """
        Cpr(T, P, state, v0, Xass0, T_step)

        Method that computes the residual heat capacity of the mixture at given
        temperature and pressure.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for the liquid phase and 'V' for the vapor phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites
        T_step: float, optional
            Step to compute the numerical temperature derivates of Helmholtz
            free energy

        Returns
        -------
        Cp: float
            residual heat capacity [J/mol K]
        """

        temp_aux = self.temperature_aux(T)
        if v0 is None:
            rho0 = None
        else:
            rho0 = 1./v0
        rho, Xass = self.density_aux(x, temp_aux, P, state, rho0, Xass0)

        rhomolecular = Na * rho

        d2a, Xass = d2ares_drho(self, x, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta

        h = T_step
        temp_aux1 = self.temperature_aux(T+h)
        temp_aux2 = self.temperature_aux(T+2*h)
        temp_aux_1 = self.temperature_aux(T-h)
        temp_aux_2 = self.temperature_aux(T-2*h)

        a1, Xass1 = dares_drho(self, x, rhomolecular, temp_aux1, Xass)
        a2, Xass2 = dares_drho(self, x, rhomolecular, temp_aux2, Xass)
        a_1, Xass_1 = dares_drho(self, x, rhomolecular, temp_aux_1, Xass)
        a_2, Xass_2 = dares_drho(self, x, rhomolecular, temp_aux_2, Xass)

        a = d2a[:2]
        da_drho = a[1] * Na
        d2a_drho = d2a[2] * Na**2

        dFdT = (a_2/12 - 2*a_1/3 + 2*a1/3 - a2/12)/h
        dFdT[1] *= Na

        d2FdT = (-a_2/12 + 4*a_1/3 - 5*a/2 + 4*a1/3 - a2/12) / h**2
        d2FdT[1] *= Na

        dP_dT = RT*(rho**2 * dFdT[1]) + P/T

        dP_drho = 2*rho*da_drho + 2.
        dP_drho += rho**2 * d2a_drho - 1.
        dP_drho *= RT

        dP_dV = -rho**2 * dP_drho
        # residual isochoric heat capacity
        Cvr_TVN = R * (-T**2*d2FdT[0] - 2*T*dFdT[0])
        # residual heat capacity
        Cpr = Cvr_TVN - R - T*dP_dT**2/dP_dV
        return Cpr

    def speed_sound(self, x, T, P, state, v0=None, Xass0=None, T_step=0.1,
                    CvId=3*R/2, CpId=5*R/2):
        """
        speed_sound(x, T, P, state, v0, Xass0, T_step)

        Method that computes the speed of sound of the mixture at a given
        temperature and pressure.

        This calculation requires that the molar weight of the fluids has been
        set in the component function.

        By default, the ideal gas Cv and Cp are set to 3R/2 and 5R/2, the user
        can supply better values if available.

        Parameters
        ----------
        x: array_like
            molar fraction array
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for the liquid phase and 'V' for the vapor phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites
        T_step: float, optional
            Step to compute the numerical temperature derivates of Helmholtz
            free energy
        CvId: float, optional
            Ideal gas isochoric heat capacity, set to 3R/2 by default [J/mol K]
        CpId: float, optional
            Ideal gas heat capacity, set to 3R/2 by default [J/mol K]

        Returns
        -------
        w: float
            speed of sound [m/s]
        """

        temp_aux = self.temperature_aux(T)
        if v0 is None:
            rho0 = None
        else:
            rho0 = 1./v0
        rho, Xass = self.density_aux(x, temp_aux, P, state, rho0, Xass0)

        rhomolecular = Na * rho

        d2a, Xass = d2ares_drho(self, x, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta

        h = T_step
        temp_aux1 = self.temperature_aux(T+h)
        temp_aux2 = self.temperature_aux(T+2*h)
        temp_aux_1 = self.temperature_aux(T-h)
        temp_aux_2 = self.temperature_aux(T-2*h)

        a1, Xass1 = dares_drho(self, x, rhomolecular, temp_aux1, Xass)
        a2, Xass2 = dares_drho(self, x, rhomolecular, temp_aux2, Xass)
        a_1, Xass_1 = dares_drho(self, x, rhomolecular, temp_aux_1, Xass)
        a_2, Xass_2 = dares_drho(self, x, rhomolecular, temp_aux_2, Xass)

        a = d2a[:2]
        da_drho = a[1] * Na
        d2a_drho = d2a[2] * Na**2

        dFdT = (a_2/12 - 2*a_1/3 + 2*a1/3 - a2/12)/h
        dFdT[1] *= Na

        d2FdT = (-a_2/12 + 4*a_1/3 - 5*a/2 + 4*a1/3 - a2/12) / h**2
        d2FdT[1] *= Na

        dP_dT = RT*(rho**2 * dFdT[1]) + P/T

        dP_drho = 2*rho*da_drho + 2.
        dP_drho += rho**2 * d2a_drho - 1.
        dP_drho *= RT

        dP_dV = -rho**2 * dP_drho
        # residual isochoric heat capacity
        Cvr_TVN = R * (-T**2*d2FdT[0] - 2*T*dFdT[0])
        # residual heat capacity
        Cpr = Cvr_TVN - R - T*dP_dT**2/dP_dV

        # speed of sound calculation
        Cp = CpId + Cpr
        Cv = CvId + Cvr_TVN

        betas = -rho * (Cv/Cp) / dP_dV

        Mwx = np.dot(x, self.Mw)
        w2 = 1000./(rho * betas * Mwx)
        w = np.sqrt(w2)

        return w
    
    def get_lnphi_pure(self, T, P, state):
        """
        get_lnphi_pure(T, P, state)
        Method that computes the logarithm of the pure component's fugacity
        coefficient at a given state, temperature T, and pressure P.
        Parameters
        ----------
        T: float
            absolute temperature [K]
        P: float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapor phase
        Returns
        -------
        lnphi_pure: float
            logarithm of pure component's fugacity coefficient
        """

        lnphi_pure = np.zeros(self.nc)
        for i, pure_eos in enumerate(self.pure_eos):
            lnphi_pure[i], _ = pure_eos.logfug(T, P, state)
        return lnphi_pure
    
    def get_lngamma(self, x, T, P, v0=None, Xass0=None, lnphi_pure=None):
        """
        get_lngamma(x, T, P, v0, Xass0)
        A method that computes the activity coefficient of the mixture at a given
        composition x, temperature T, and pressure P.
        Parameters
        ----------
        x: array_like
            molar fraction array
        T: float
            absolute temperature [K]
        P: float
            pressure [Pa]
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of the fraction of non-bonded sites
        lnphi_pure: array, optional
            logarithm of the pure components's fugacity coefficient.
            Computed if not provided.
        Returns
        -------
        lngamma: float
            logarithm of activity coefficient model
        """


        if isinstance(x, (float, int)):
            x = np.array([x, 1.-x])
        elif not isinstance(x, np.ndarray):
            x = np.array(x)

        if self.nc > 2 and x.shape[0] < 2:
            raise ValueError('Please supply the whole molfrac vector for non-binary mixtures')

        lnphi_mix, _ = self.logfugef(x, T, P, 'L', v0=v0, Xass0=Xass0)

        if lnphi_pure is None:
            lnphi_pure = np.zeros_like(x)
            for i, pure_eos in enumerate(self.pure_eos):
                lnphi_pure[i], _ = pure_eos.logfug(T, P, 'L')

        lngamma = lnphi_mix - lnphi_pure
        return lngamma
    
    def contributions(self, x, rho, T, Xass0=None):
        temp_aux = self.temperature_aux(T)
        out = all_contributions(self, x, rho, temp_aux, Xass0)
        Mono, Chain, Disp,  Asso, Polar, DH, Born = out
        return Mono, Chain, Disp,  Asso, Polar, DH, Born

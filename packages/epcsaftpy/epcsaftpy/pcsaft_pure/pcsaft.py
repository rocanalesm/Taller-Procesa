from __future__ import division, print_function, absolute_import
import numpy as np
from .ideal.ideal import aideal, daideal_drho, d2aideal_drho
from .association.association_aux import association_config
from .dispersion.coefficients import aji, bji 

from .ares import ares, dares_drho, d2ares_drho, all_contributions


from .routines.density_solver import density_newton
from .routines.psat_saft import psat
from .routines.tsat_saft import tsat
from .routines.critical_pure import get_critical

from ..constants import kb, Na


R = Na * kb




class pcsaft_pure():
    '''
    Pure component PC-SAFT EoS Object
    
    This object have implemeted methods for phase equilibrium
    as for interfacial properties calculations.

    Parameters
    ----------
    pure : object
        pure component created with component class
    compute_critical: bool
        If True the critical point of the fluid will attempt to be computed
        (it might fail for some fluids).

    Attributes
    ----------
    ms: number of chain segments
    sigma: segment diameter parameter [m]
    eps: dispersion energy [J]

    eABij: association energy [J]
    kappaABij: association volume
    sites: triplet of number of association sites [B, P, N]

    mupol: dipolar moment [Debye]
    xpol: fraction of dipolar segment on a chain molecule

    cii : influence factor for SGT [J m^5 / mol^2]

    Methods
    -------
    temperature_aux : computes temperature dependent parameters of the fluid
    density : computes the density of the fluid
    psat : computes saturation pressure
    tsat : computes saturation temperature
    get_critical : attemps to compute the critical point of the fluid
    afcn: computes total Helmholtz energy
    dafcn_drho : computes total Helmholtz energy and its density derivative
    d2afcn_drho : computes total Helmholtz energy and it density derivatives
    pressure : computes the pressure
    dP_drho : computes pressure and its density derivative
    logfug : computes the fugacity coefficient
    a0ad : computes adimentional Helmholtz density energy
    muad : computes adimentional chemical potential
    dOm : computes adimentional Thermodynamic Grand Potential
    ci :  computes influence parameters matrix for SGT
    sgt_adim : computes adimentional factors for SGT

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
    logfug_aux : computes logfug
    a0ad_aux : compute a0ad
    muad_aux : computes muad
    dOm_aux : computes dOm
    
    
    Code template based in SGTpy (https://github.com/gustavochm/sgtpy)
    '''

    def __init__(self, pure, compute_critical=True):

        self.pure = pure
        self.Mw = pure.Mw
        self.ms = pure.ms
        self.sigmaT_bool = pure.sigmaT_bool
        if not self.sigmaT_bool:
            self.sigma = pure.sigma0
            self.sigma3 = self.sigma**3
        else:
            self.sigma = pure.sigma0
            self.sigma0 = pure.sigma0
            self.t1 = pure.t1
            self.t2 = pure.t2
            self.t3 = pure.t3
            self.t4 = pure.t4
        self.eps = pure.eps

        # coefficients
        mm1 = (self.ms - 1)/self.ms
        mm2 = mm1*(self.ms - 2)/self.ms
        self.ai = aji[0] + aji[1]*mm1 + aji[2]*mm2
        self.bi = bji[0] + bji[1]*mm1 + bji[2]*mm2
       

        # association pre-configuration
        self.eABij = pure.eAB
        self.kappaABij = pure.kappaAB
        

        self.sites = pure.sites
        S, DIJ, indexabij, nsites, diagasso = association_config(self)
        assoc_bool = self.eABij != 0
        self.assoc_bool = assoc_bool
        if assoc_bool:
            self.S = S
            self.DIJ = DIJ
            self.indexabij = indexabij
            self.nsites = nsites
            self.diagasso = diagasso

        # polar pre-configuration
        self.mupol = pure.mupol
        self.xpol = pure.xpol
        polar_bool = self.xpol != 0
        self.polar_bool = polar_bool
        if polar_bool:
            # 1 D = 3.33564e-30 C * m
            # 1 C^2 = 9e9 N m^2
            cte = (3.33564e-30)**2 * (9e9)
            self.xpmupol2 = self.xpol*self.mupol**2*cte
        else:
            self.xpmupol2 = 0


        # for SGT computations
        self.cii = np.array(pure.cii, ndmin=1)

        # fir entropy scaling calculation
        self.viscosity_parameters = np.array(pure.viscosity_parameters)

        # computing critical point
        self.critical = False
        if compute_critical:
            try:
                out = get_critical(self, None, None, method='hybr',
                                full_output=True)
                if out.success:
                    self.critical = True
                    self.Tc = out.Tc
                    self.Pc = out.Pc
                    self.rhoc = out.rhoc
            except:
                print('The critical point could not be computed.')




    
    def temperature_aux(self, T):
        """
        temperature_aux(T)

        Method that computes temperature dependent parameters.
        It returns the following list:

        temp_aux = [beta, eps_beta, d, mes3, m2e2s3, aux_xi]



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
            self.sigma3 = self.sigma**3
        beta = 1 / (kb*T)
        eps_beta = self.eps*beta
        d = self.sigma * (1 - 0.12 * np.exp(-3 * eps_beta))
        dia3 = d**3

        # Parameters needed for evaluating the helmothlz contributions

        mes3 = self.ms * self.sigma**3 * eps_beta
        m2e2s3 = mes3 * self.ms * eps_beta

        deta_drho = self.ms * np.pi * dia3 / 6
        Fab = np.exp(beta * self.eABij) - 1.
        
        betaxpmu2 = beta*self.xpmupol2
        
        temp_aux = [beta, eps_beta, d, dia3, mes3, m2e2s3, deta_drho, Fab, betaxpmu2]
        return temp_aux

    def density_aux(self, temp_aux, P, state, rho0=None, Xass0=None):
        """
        density_aux(T, temp_aux, state, rho0, Xass0)
        Method that computes the density of the fluid at T, P

        Parameters
        ----------
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapor phase
        rho0 : float, optional
            initial guess to compute density root [mol/m^3]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        density: float
            density [mol/m^3]
        Xass : array
            computed fraction of nonbonded sites

        """
        etamax = 0.7404804896930609
        rhomax = (6 * etamax) / (self.ms * np.pi * self.sigma**3) / Na
        if rho0 is None:
            if state == 'L':
                eta0 = 0.5 * P**0
                rho0 = rhomax * eta0 / etamax
            elif state == 'V':
                #eta0 = 1e-10
                #rho0 = (6 * eta0) / (self.ms * np.pi * self.sigma**3) / Na
                beta = temp_aux[0]
                rho0 = P * beta / Na


        rho, Xass = density_newton(rho0, temp_aux, P, Xass0, rhomax, state, self)
        return rho, Xass

    def density(self, T, P, state, rho0=None, Xass0=None):
        """
        density(T, P, state)
        Method that computes the density of the fluid at T, P

        Parameters
        ----------

        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapor phase
        rho0 : float, optional
            initial guess to compute density root [mol/m^3]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        density: float
            density [mol/m^3]
        """
        temp_aux = self.temperature_aux(T)
        rho, Xass = self.density_aux(temp_aux, P, state, rho0, Xass0)
        return float(rho)

    def psat(self, T, P0=None, v0=[None, None], Xass0=[None, None],
             full_output=False):
        """
        psat(T, P0)

        Method that computes saturation pressure at fixed T

        Parameters
        ----------

        T : float
            absolute temperature [K]
        P0 : float, optional
            initial value to find saturation pressure [Pa]
        v0: list, optional
            initial guess for liquid and vapor phase, respectively [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites
        full_output: bool, optional
            whether to outputs or not all the calculation info.

        Returns
        -------
        psat : float
            saturation pressure [Pa]
        vl : float
            liquid saturation volume [m3/mol]
        vv : float
            vapor saturation volume [m3/mol]
        """
        out = psat(self, T, P0, v0, Xass0, full_output)
        return out

    def tsat(self, P,  T0=None, Tbounds=None, v0=[None, None],
             Xass0=[None, None], full_output=False):
        """
        tsat(P, Tbounds)

        Method that computes saturation temperature at given pressure.

        Parameters
        ----------

        P : float
            absolute pressure [Pa]
        T0 : float, optional
             Temperature to start iterations [K]
        Tbounds : tuple, optional
                (Tmin, Tmax) Temperature interval to start iterations [K]
        v0: list, optional
            initial guess for liquid and vapor phase, respectively [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites
        full_output: bool, optional
            whether to outputs or not all the calculation info.

        Returns
        -------
        tsat : float
            saturation temperature [K]
        vl : float
            liquid saturation volume [m^3/mol]
        vv : float
            vapor saturation volume [m^3/mol]
        """
        out = tsat(self, P, T0, Tbounds, v0, Xass0, full_output)
        return out

    def get_critical(self, Tc0=None, rhoc0=None, method='hybr',
                     full_output=False, overwrite=False):
        """
        get_critical(Tc0, rhoc0, method)

        Method that solves the critical coordinate of the fluid.
        This metho requires good initial guesses for the critical temperature
        and density to converge.

        Second derivative of pressure against volume is estimated numerically.

        Parameters
        ----------
        Tc0 : float, optional
            initial guess for critical temperature [K]
        rhoc : float, optional
            initial guess for critical density [mol/m^3]
        method : string, optional
            SciPy; root method to solve critical coordinate
        full_output: bool, optional
            whether to outputs or not all the calculation info
        overwrite: bool, optional
            wheter to overwrite already computed critical points

        Returns
        -------
        Tc: float
            Critical temperature [K]
        Pc: float
            Critical pressure [Pa]
        rhoc: float
            Critical density [mol/m3]
        """
        out = get_critical(self, Tc0, rhoc0, method, full_output)
        if overwrite:
            if full_output:
                if out.success:
                    self.critical = True
                    self.Tc = out.Tc
                    self.Pc = out.Pc
                    self.rhoc = out.rhoc
            else:
                Tc0 = out[0]
                rhoc0 = out[2]
                out2 = get_critical(self, Tc0, rhoc0, method, full_output=True)
                if out2.success:
                    self.critical = True
                    self.Tc = out2.Tc
                    self.Pc = out2.Pc
                    self.rhoc = out2.rhoc
        return out

    def ares(self, rho, T, Xass0=None):
        """
        ares(x, rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the fluid.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: float
           residual dimentionless Helmholtz free energy [Adim]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = ares(self, rho, temp_aux, Xass0)
        return a, Xass

    def dares_drho(self, rho, T, Xass0=None):
        """
        dares_drho(rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the fluid
        and its first density derivative.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array_like
           residual dimentionless Helmholtz free energy [Adim, m^3]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = dares_drho(self, rho, temp_aux, Xass0)
        return a, Xass

    def d2ares_drho(self, rho, T, Xass0=None):
        """
        d2ares_drho(rho, T, Xass0)
        Method that computes the residual Helmholtz free energy of the fluid
        and its first and second density derivatives.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T: float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array_like
           residual dimentionless Helmholtz free energy [Adim, m^3, m^6]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = d2ares_drho(self, rho, temp_aux, Xass0)
        return a, Xass

    def afcn_aux(self, rho, temp_aux, Xass0=None):
        """
        afcn_aux(rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the fluid.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: float
           Helmholtz free energy [J/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        a, Xass = ares(self, rho, temp_aux, Xass0)
        a += aideal(rho, beta)
        a *= (Na/beta)
        return a, Xass

    def dafcn_aux(self, rho, temp_aux, Xass0=None):
        """
        dafcn_aux(rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the fluid and
        its first density derivative.

        Parameters
        ----------
        rho: float
            density [mol/m3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array
           Helmholtz free energy and its derivative  [J/mol, J m^3/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        a, Xass = dares_drho(self, rho, temp_aux, Xass0)
        a += daideal_drho(rho, beta)
        a *= (Na/beta)
        return a, Xass

    def d2afcn_aux(self, rho, temp_aux, Xass0=None):
        """
        d2afcn_aux(rho, temp_aux, Xass0)
        Method that computes the total Helmholtz free energy of the fluid and
        its first ans second density derivative.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array
           Helmholtz free energy and its derivatives: a, da, d2a
           [J/mol, J m^3/mol^2,  J m^6/mol^3]
        Xass : array
            computed fraction of nonbonded sites
        """
        beta = temp_aux[0]
        a, Xass = d2ares_drho(self, rho, temp_aux, Xass0)
        a += d2aideal_drho(rho, beta)
        a *= (Na/beta)
        return a, Xass

    def afcn(self, rho, T, Xass0=None):
        """
        afcn(rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the fluid.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: float
           Helmholtz free energy [J/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = self.afcn_aux(rho, temp_aux, Xass0)
        return a

    def dafcn_drho(self, rho, T, Xass0=None):
        """
        dafcn_drho(rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the fluid and
        its first density derivative.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array
           Helmholtz free energy and its derivative  [J/mol, J m^3/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = self.dafcn_aux(rho, temp_aux, Xass0)
        return a

    def d2afcn_drho(self, rho, T, Xass0=None):
        """
        d2afcn_drho(rho, T, Xass0)
        Method that computes the total Helmholtz free energy of the fluid and
        its first ans second density derivative.

        Parameters
        ----------
        rho: float
            molecular density [molecules/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a: array
           Helmholtz free energy and its derivatives: a, da, d2a
           [J/mol, J m^3/mol,  J m^6/mol]
        """
        temp_aux = self.temperature_aux(T)
        a, Xass = self.d2afcn_aux(rho, temp_aux, Xass0)
        return a

    def pressure_aux(self, rho, temp_aux, Xass0=None):
        """
        pressure_aux(rho, temp_aux, Xass0)

        Method that computes the pressure at given density [mol/m3] and
        temperature [K]

        Parameters
        ----------
        rho: float
            density [mol/m3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        Xass : array
            computed fraction of nonbonded sites
        """
        rhomolecular = Na * rho
        da, Xass = self.dafcn_aux(rhomolecular, temp_aux, Xass0)
        afcn, dafcn = da
        Psaft = rhomolecular**2 * dafcn / Na
        return Psaft, Xass

    def dP_drho_aux(self, rho, temp_aux, Xass0=None):
        """
        dP_drho_aux(rho, temp_aux, Xass0)

        Method that computes the pressure and its density derivative at given
        density [mol/m3] and temperature [K]

        Parameters
        ----------
        rho: float
            density [mol/m3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

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
        da, Xass = self.d2afcn_aux(rhomolecular, temp_aux, Xass0)
        afcn, dafcn, d2afcn = da
        Psaft = rhomolecular**2 * dafcn / Na
        dPsaft = 2 * rhomolecular * dafcn + rhomolecular**2 * d2afcn
        return Psaft, dPsaft, Xass

    def pressure(self, rho, T, Xass0=None):
        """
        pressure(rho, T, Xass0)

        Method that computes the pressure at given density [mol/m3] and
        temperature [K]

        Parameters
        ----------
        rho: float
            density [mol/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        """
        temp_aux = self.temperature_aux(T)
        Psaft, Xass = self.pressure_aux(rho, temp_aux, Xass0)
        return Psaft

    def dP_drho(self, rho, T, Xass0=None):
        """
        dP_drho(rho, T, Xass0)

        Method that computes the pressure and its density derivative at given
        density [mol/m3] and temperature [K]

        Parameters
        ----------
        rho: float
            density [mol/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        P : float
            pressure [Pa]
        dP: float
            derivate of pressure respect density [Pa m^3 / mol]
        """
        temp_aux = self.temperature_aux(T)
        Psaft, dPsaft, Xass = self.dP_drho_aux(rho, temp_aux, Xass0)

        return Psaft, dPsaft

    def logfug_aux(self, temp_aux, P, state, v0=None, Xass0=None):
        """
        logfug_aux(T, P, state, v0, Xass0)

        Method that computes the fugacity coefficient at given
        composition, temperature and pressure.

        Parameters
        ----------
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        logfug : float
            fugacity coefficient
        v : float
            computed volume of the phase [m^3/mol]
        Xass : array
            computed fraction of nonbonded sites
        """
        if v0 is None:
            rho, Xass = self.density_aux(temp_aux, P, state, None, Xass0)
        else:
            rho0 = 1./v0
            rho, Xass = self.density_aux(temp_aux, P, state, rho0, Xass0)
        v = 1./rho
        rhomolecular = Na * rho
        ar, Xass = ares(self, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta
        Z = P * v / RT
        lnphi = ar + (Z - 1.) - np.log(Z)
        return lnphi, v, Xass

    def logfug(self, T, P, state, v0=None, Xass0=None):
        """
        logfug(T, P, state, v0, Xass0)

        Method that computes the fugacity coefficient at given temperature
        and pressure.

        Parameters
        ----------
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        logfug: float
            fugacity coefficient
        v: float
            computed volume of the phase [m^3/mol]
        """
        temp_aux = self.temperature_aux(T)
        lnphi, v, Xass = self.logfug_aux(temp_aux, P, state, v0, Xass0)
        return lnphi, v

    def ci(self, T):
        '''
        ci(T)

        Method that evaluates the polynomial for the influence parameters used
        in the SGT theory for surface tension calculations.

        Parameters
        ----------
        T : float
            absolute temperature [K]

        Returns
        -------
        ci: float
            influence parameters [J m5 mol-2]
        '''

        return np.polyval(self.cii, T)

    def sgt_adim(self, T):
        '''
        sgt_adim(T)

        Method that evaluates adimentional factor for temperature, pressure,
        density, tension and distance for interfacial properties computations
        with SGT.

        Parameters
        ----------
        T : float
        absolute temperature [K]

        Returns
        -------
        Tfactor : float
            factor to obtain dimentionless temperature (K -> K)
        Pfactor : float
            factor to obtain dimentionless pressure (Pa -> Pa/RT)
        rofactor : float
            factor to obtain dimentionless density (mol/m3 -> mol/m3)
        tenfactor : float
            factor to obtain dimentionless surface tension (mN/m)
        zfactor : float
            factor to obtain dimentionless distance  (Amstrong -> m)
        '''
        cii = self.ci(T)  # computing temperature dependent cii

        Tfactor = 1.
        Pfactor = 1.
        rofactor = 1.
        tenfactor = np.sqrt(cii) * 1000  # To give tension in mN/m
        zfactor = 10**-10

        return Tfactor, Pfactor, rofactor, tenfactor, zfactor

    def sgt_adim_fit(self, T):

        Tfactor = 1
        Pfactor = 1
        rofactor = 1
        tenfactor = 1. * 1000  # To give tension in mN/m

        return Tfactor, Pfactor, rofactor, tenfactor

    def a0ad_aux(self, rho, temp_aux, Xass0=None):
        """
        a0ad_aux(ro, temp_aux, Xass0)

        Method that computes the adimenstional Helmholtz density energy at
        given density and temperature.

        Parameters
        ----------

        rho : float
            density [mol/m^3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a0ad: float
            Helmholtz density energy [J/m^3]
        Xass : array
            computed fraction of nonbonded sites
        """
        rhomolecular = rho * Na
        a0, Xass = self.afcn_aux(rhomolecular, temp_aux, Xass0)
        a0 *= rho

        return a0, Xass

    def a0ad(self, rho, T, Xass0=None):
        """
        a0ad(ro, T, Xass0)

        Method that computes the adimenstional Helmholtz density energy at
        given density and temperature.

        Parameters
        ----------

        rho : float
            density [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        a0ad: float
            Helmholtz density energy [J/m^3]
        """
        temp_aux = self.temperature_aux(T)
        a0, Xass = self.a0ad_aux(rho, temp_aux, Xass0)
        return a0

    def muad_aux(self, rho, temp_aux, Xass0=None):
        """
        muad_aux(rho, temp_aux, Xass0)

        Method that computes the adimenstional chemical potential at given
        density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        muad: float
            chemical potential [J/mol]
        Xass : array
            computed fraction of nonbonded sites
        """

        rhomolecular = rho * Na
        da, Xass = self.dafcn_aux(rhomolecular, temp_aux, Xass0)
        afcn, dafcn = da
        mu = afcn + rhomolecular * dafcn

        return mu, Xass

    def muad(self, rho, T, Xass0=None):
        """
        muad(rho, T, Xass0)

        Method that computes the adimenstional chemical potential at given
        density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        muad: float
            chemical potential [J/mol]
        """
        temp_aux = self.temperature_aux(T)
        mu, Xass = self.muad_aux(rho, temp_aux, Xass0)
        return mu

    def dOm_aux(self, rho, temp_aux, mu, Psat, Xass0=None):
        """
        dOm_aux(rho, temp_aux, mu, Psat, Xass0)

        Method that computes the adimenstional Thermodynamic Grand potential
        at given density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        temp_aux : list
            temperature dependend parameters computed with temperature_aux(T)
        mu : float
            adimentional chemical potential at equilibrium
        Psat : float
            adimentional pressure [Pa]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        GPT: float
            Thermodynamic Grand potential [Pa]
        Xass : array
            computed fraction of nonbonded sites
        """
        a0, Xass = self.a0ad_aux(rho, temp_aux, Xass0)
        GPT = a0 - rho*mu + Psat

        return GPT, Xass

    def dOm(self, rho, T, mu, Psat, Xass0=None):
        """
        dOm(rho, T, mu, Psat, Xass0)

        Method that computes the adimenstional Thermodynamic Grand potential
        at given density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        T : float
            absolute temperature [K]
        mu : float
            adimentional chemical potential at equilibrium
        Psat : float
            adimentional pressure [Pa]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        Out: float
            Thermodynamic Grand potential [Pa]
        """
        temp_aux = self.temperature_aux(T)
        GPT, Xass = self.dOm_aux(rho, temp_aux, mu, Psat, Xass0)
        return GPT

    def EntropyR(self, T, P, state, v0=None, Xass0=None):
        """
        EntropyR(T, P, state, v0, Xass0)

        Method that computes the residual entropy (NPT) at given temperature and
        pressure.

        Parameters
        ----------
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        Sr : float
            residual entropy [J/mol K]

        """
        temp_aux = self.temperature_aux(T)
        if v0 is None:
            rho, Xass = self.density_aux(temp_aux, P, state, None, Xass0)
        else:
            rho0 = 1./v0
            rho, Xass = self.density_aux(temp_aux, P, state, rho0, Xass0)
        
        v = 1./rho
        beta = temp_aux[0]
        RT = Na/beta
        rhomolecular = Na * rho
        Z = P * v / RT  

        a, Xass = ares(self, rhomolecular, temp_aux, Xass0)

        if type(Xass) == np.ndarray:
            Xass = Xass.astype('complex128')

        h = np.finfo(1.0).eps
        temp_aux_h = self.temperature_aux(T + h * 1j )

        ah, _ = ares(self, rhomolecular + 0. * 1j, temp_aux_h, Xass)


        dadT = ah.imag/h
        Sr_TVN = -T*dadT - a  # residual entropy (TVN) divided by R
        Sr_TPN = Sr_TVN + np.log(Z)  # residual entropy (TPN) divided by R
        Sr_TPN *= R
        return Sr_TPN

    def EntropyR_NrhoT(self, rho, T, Xass0 = None):
        """
        EntropyR_NrhoT(rho, T, Xass0 = None)

        Method that computes the residual entropy (NVT) at given density and
        temperature.

        Parameters
        ----------
        rho: float
            density [mol/m3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites


        Returns
        -------
        Sr : float
            residual entropy [J/mol K]

        """
        temp_aux = self.temperature_aux(T)
        rhomolecular = Na * rho

        a, Xass = ares(self, rhomolecular, temp_aux, Xass0)
        if type(Xass) == np.ndarray:
            Xass = Xass.astype('complex128')

        h = np.finfo(1.0).eps
        temp_aux_h = self.temperature_aux(T + h * 1j )

        ah, _ = ares(self, rhomolecular + 0. * 1j, temp_aux_h, Xass)


        dadT = ah.imag/h
        Sr_TVN = -T * dadT - a  # residual entropy (TVN) divided by R
        Sr_TVN *= R
        return Sr_TVN

    def EnthalpyR(self, T, P, state, v0=None, Xass0=None):
        """
        EnthalpyR(T, P, state, v0, Xass0, T_step)

        Method that computes the residual enthalpy at given temperature and
        pressure.

        Parameters
        ----------
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites

        Returns
        -------
        Hr : float
            residual enthalpy [J/mol]

        """
        temp_aux = self.temperature_aux(T)
        if v0 is None:
            rho, Xass = self.density_aux(temp_aux, P, state, None, Xass0)
        else:
            rho0 = 1./v0
            rho, Xass = self.density_aux(temp_aux, P, state, rho0, Xass0)
        v = 1./rho
        rhomolecular = Na * rho
        beta = temp_aux[0]
        RT = Na/beta
        Z = P * v / RT  

        a, Xass = ares(self, rhomolecular, temp_aux, Xass)

        if type(Xass) == np.ndarray:
            Xass = Xass.astype('complex128')
            
        h = np.finfo(1.0).eps
        temp_aux_h = self.temperature_aux(T + h * 1j )

        ah, _ = ares(self, rhomolecular + 0. * 1j, temp_aux_h, Xass)

        dadT = ah.imag/h
        Sr_TVN = -T*dadT - a  # residual entropy divided by R
        Hr_TPN = a + Sr_TVN + Z - 1.  # residual entalphy divided by RT
        Hr_TPN *= RT
        return Hr_TPN

    def CvR(self, rho, T, Xass0=None, T_step=0.1):
        """
        CvR(rho, T, Xass0, T_step)

        Method that computes the residual isochoric heat capacity at given
        density and temperature.

        Parameters
        ----------
        rho : float
            density [mol/m^3]
        T : float
            absolute temperature [K]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites
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

        a, Xass = ares(self, rhomolecular, temp_aux, Xass0)

        h = T_step
        temp_aux1 = self.temperature_aux(T+h)
        temp_aux2 = self.temperature_aux(T+2*h)
        temp_aux_1 = self.temperature_aux(T-h)
        temp_aux_2 = self.temperature_aux(T-2*h)

        a1, Xass1 = ares(self, rhomolecular, temp_aux1, Xass)
        a2, Xass2 = ares(self, rhomolecular, temp_aux2, Xass)
        a_1, Xass_1 = ares(self, rhomolecular, temp_aux_1, Xass)
        a_2, Xass_2 = ares(self, rhomolecular, temp_aux_2, Xass)

        dFdT = (a_2/12 - 2*a_1/3 + 2*a1/3 - a2/12)/h
        d2FdT = (-a_2/12 + 4*a_1/3 - 5*a/2 + 4*a1/3 - a2/12)/h**2

        Cvr_TVN = -T**2*d2FdT - 2*T*dFdT  # residual isochoric heat capacity
        Cvr_TVN *= R
        return Cvr_TVN

    def CpR(self, T, P, state, v0=None, Xass0=None, T_step=0.1):
        """
        Cpr(T, P, state, v0, Xass0, T_step)

        Method that computes the residual heat capacity at given temperature
        and pressure.

        Parameters
        ----------
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites
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
            rho, Xass = self.density_aux(temp_aux, P, state, None, Xass0)
        else:
            rho0 = 1./v0
            rho, Xass = self.density_aux(temp_aux, P, state, rho0, Xass0)

        rhomolecular = Na * rho

        d2a, Xass = d2ares_drho(self, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta

        h = T_step
        temp_aux1 = self.temperature_aux(T+h)
        temp_aux2 = self.temperature_aux(T+2*h)
        temp_aux_1 = self.temperature_aux(T-h)
        temp_aux_2 = self.temperature_aux(T-2*h)

        a1, Xass1 = dares_drho(self, rhomolecular, temp_aux1, Xass)
        a2, Xass2 = dares_drho(self, rhomolecular, temp_aux2, Xass)
        a_1, Xass_1 = dares_drho(self, rhomolecular, temp_aux_1, Xass)
        a_2, Xass_2 = dares_drho(self, rhomolecular, temp_aux_2, Xass)

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

    def speed_sound(self, T, P, state, v0=None, Xass0=None, T_step=0.1,
                    CvId=3*R/2, CpId=5*R/2):
        """
        speed_sound(T, P, state, v0, Xass0, T_step, CvId, CpId)

        Method that computes the speed of sound at given temperature
        and pressure.

        This calculation requires that the molar weight of the fluid has been
        set in the component function.

        By default the ideal gas Cv and Cp are set to 3R/2 and 5R/2, the user
        can supply better values if available.

        Parameters
        ----------
        T : float
            absolute temperature [K]
        P : float
            pressure [Pa]
        state : string
            'L' for liquid phase and 'V' for vapour phase
        v0: float, optional
            initial guess for volume root [m^3/mol]
        Xass0: array, optional
            Initial guess for the calculation of fraction of non-bonded sites
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
            rho, Xass = self.density_aux(temp_aux, P, state, None, Xass0)
        else:
            rho0 = 1./v0
            rho, Xass = self.density_aux(temp_aux, P, state, rho0, Xass0)

        rhomolecular = Na * rho

        d2a, Xass = d2ares_drho(self, rhomolecular, temp_aux, Xass)
        beta = temp_aux[0]
        RT = Na/beta

        h = T_step
        temp_aux1 = self.temperature_aux(T+h)
        temp_aux2 = self.temperature_aux(T+2*h)
        temp_aux_1 = self.temperature_aux(T-h)
        temp_aux_2 = self.temperature_aux(T-2*h)

        a1, Xass1 = dares_drho(self, rhomolecular, temp_aux1, Xass)
        a2, Xass2 = dares_drho(self, rhomolecular, temp_aux2, Xass)
        a_1, Xass_1 = dares_drho(self, rhomolecular, temp_aux_1, Xass)
        a_2, Xass_2 = dares_drho(self, rhomolecular, temp_aux_2, Xass)

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

        w2 = 1000./(rho * betas * self.Mw)
        w = np.sqrt(w2)

        return w
    
    def contributions(self, rho, T, Xass0=None):
        temp_aux = self.temperature_aux(T)
        Mono, Chain, Disp,  Asso, Polar = all_contributions(self, rho, temp_aux, Xass0)
        return Mono, Chain, Disp,  Asso, Polar

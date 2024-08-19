from __future__ import division, print_function, absolute_import
import numpy as np
from copy import copy
from .constants import kb
import json
import pandas as pd
from IPython.display import HTML, display


class component(object):
    '''
    Creates an object with pure component info

    Parameters
    ----------
    name : str
        Name of the component
    ms : float
        Chain lenght
    sigma : float or list
        Segment diameter, input in [Ãngstrom], stored in [m]
    eps : float
        Dispersion energy, input in [k], stored in [J]
    eAB : float
        Association Energy, input in [k], stores in [J]
    kappaAB : float
        Association volume
    sites : list
        Association sites [Bipolar, Positive, Negative]
    mupol : float
        dipolar moment [Debye]
    xpol : float
        fraction of dipolar segment on a chain molecule
        (Jog and Chapman polar term)
    z : float
        ions valence
    er : float
        dielectric constant
    Mw : float
        molar weight [g/mol]
    cii : list
        polynomial coefficient for influence parameter
        in SGT for sgtpy module [J m^5 / mol^2]
    '''

    def __init__(self, name='None', ms=1., sigma=0., eps=0., eAB=0., 
                 kappaAB=0, sites=[0, 0, 0], mupol=0, xpol=0., z = 0, er = 0, 
                 Mw=1., cii=0., viscosity_parameters = [0, 0, 0, 0], 
                 pure_path = None, cas=None, smiles=None, option_name='name'):

        self.name = name
        self.cii = np.atleast_1d(cii)  # Influence factor SGT, list or array
        self.Mw = Mw  # molar weight in g/mol
        self.nc = 1

        # SAFT parameters
        self.ms = ms
        self.sigmaT_bool = type(sigma)==list
        
        if self.sigmaT_bool:
            self.sigma0 = sigma[0] * 1e-10  # meters
            self.t1 = sigma[1] * 1e-10      # meters
            self.t2 = sigma[2]              # 1/K
            self.t3 = sigma[3] * 1e-10      # meters
            self.t4 = sigma[4]              # 1/K
        else:
            self.sigma0 = sigma * 1e-10  # meters
            self.t1 = 0
            self.t2 = 0
            self.t3 = 0
            self.t4 = 0
        
        self.eps = eps * kb  # Joule

        # Association parameters
        self.eAB = eAB * kb  # Joule
        self.kappaAB = kappaAB 
        self.sites = sites

        # Polar parameters
        self.mupol = mupol  # Debye
        self.xpol = xpol
        
        # Electrolyte parameters 
        self.z = z
        self.erT_bool = type(er)==list
        if self.erT_bool: 
            # er = er0 + er1*T + er2*ln T
            self.er0 = er[0]
            self.er1 = er[1]
            self.er2 = er[2]
        else:
            self.er0 = er
            self.er1 = 0.
            self.er2 = 0.
        self.viscosity_parameters = viscosity_parameters

        self.reference = 'This work'

        if pure_path != None:
            f = open(pure_path)
            data = json.load(f)
            Dcas = []
            Dname = []
            Diupac = []
            Dsimles = []
            Dinchi = []
            Dformula = []
            for i in range(len(data)):
                dat = data[i]["identifier"]
                Dcas.append(dat["cas"])
                Dname.append(dat["name"])
                try:
                    Diupac.append(dat["iupac_name"])
                except:
                    Diupac.append("None")
                Dsimles.append(dat["smiles"])
                Dinchi.append(dat["inchi"])
                try:
                    Dformula.append(dat["formula"])
                except:
                    Dformula.append("None")

            if option_name=="name":
                index = Dname.index(name)
            elif option_name=="cas":
                index = Dcas.index(cas)
            elif option_name=="smiles":
                index = Dsimles.index(smiles)
            else:
                raise Exception("Warning. That 'option_name' is not implemented.")
            dat = data[index]["identifier"]
            self.cas = dat["cas"]
            self.smiles = dat["smiles"]
            self.Mw = data[index]["molarweight"]

            dat_parameters = data[index]["model_record"]
            self.ms = float(dat_parameters["m"])
            self.sigma0 = float(dat_parameters["sigma"])*1e-10
            try:
                self.sigmaT_bool = True
                self.t1 = float(dat_parameters["t1"]) * 1e-10      # meters
                self.t2 = float(dat_parameters["t2"]) 
                self.t3 = float(dat_parameters["t3"]) * 1e-10      # meters
                self.t4 = float(dat_parameters["t4"])
            except:
                self.sigmaT_bool = False
                self.t1 = 0
                self.t2 = 0
                self.t3 = 0
                self.t4 = 0

            self.eps = float(dat_parameters["epsilon_k"]) * kb  

            # Association parameters
            try:
                self.eAB = float(dat_parameters["epsilon_k_ab"]) * kb  
                self.kappaAB = float(dat_parameters["kappa_ab"])
                try:
                    na = int(dat_parameters["na"]) 
                    nb = int(dat_parameters["nb"]) 
                    try:
                        nc = int(dat_parameters["nc"]) # Bivalent
                    except:
                        nc = 0
                    self.sites = [nc, na, nb]
                except:
                    self.sites = [0, 1, 1]             # 2B
            except:
                self.eAB = 0
                self.kappaAB = 0  

            # Polar parameters
            try:
                self.mupol = float(dat_parameters["mu"])
            except:
                pass      

            # Electrolyte parameters
            try:
                er = dat_parameters["er"]
                self.erT_bool = type(er)==list
                if self.erT_bool: 
                    # er = er0 + er1*T + er2*ln T
                    self.er0 = er[0]
                    self.er1 = er[1]
                    self.er2 = er[2]
                else:
                    self.er0 = er
                    self.er1 = 0.
                    self.er2 = 0.
            except:
                pass
            
            try:
                self.z = float(dat_parameters["z"])
            except:
                pass

            # Viscosity parameters
            try:
                self.viscosity_parameters = dat_parameters["viscosity"]
            except:
                pass   

            index = pure_path.find('/')
            n = len(pure_path)
            self.reference = pure_path[index+1:n]
            


    def __add__(self, component2):
        '''
        Methods to add two components and create a mixture with them.
        '''
        return mixture(self, component2)

    def ci(self, T):
        """
        Method that evaluates the polynomial for cii coeffient of SGT for
        sgtpy module 
        cii must be in J m^5 / mol^2 and T in K.

        Parameters
        ----------
        T : float
            absolute temperature [K]

        Returns
        -------
        ci : float
            influence parameter at given temperature [J m^5 / mol^2]

        """

        return np.polyval(self.cii, T)
    
    def printParameters(self):
        """
        Method that prints the PC-SAFT parameters for a component.
        """
        pcsaft = {'component': [self.name], 'Mw': [self.Mw], '$ms$': [self.ms], '$\sigma$': [self.sigma0*1e10], r'$\epsilon / k_B$': [self.eps/kb]}
        if self.sigmaT_bool:
            pcsaft = {'component': [self.name], 'Mw': [self.Mw], '$ms$': [self.ms], '$\sigma$': ['$\sigma(T)$']}
        df_pcsaft = pd.DataFrame(data=pcsaft)
        if max(self.sites) != 0:
            assoc = {'[B, P, N]': [self.sites], r'$\epsilon _{AB} / k_B$': [self.eAB/kb], r'$k_{AB}$': [self.kappaAB]}
            df_assoc = pd.DataFrame(data=assoc)
            df_pcsaft = pd.concat([df_pcsaft, df_assoc], axis=1)
        if self.xpol != 0:
            polar = {'$x_p$': [self.xpol], r'$\mu$': [self.mupol]}
            df_polar = pd.DataFrame(data=polar)
            df_pcsaft = pd.concat([df_pcsaft, df_polar], axis=1)
        if self.z != 0:
            if self.erT_bool:
                electrolyte = {'$z$': [self.z], r'$\epsilon_r$': ['$\epsilon_r(T)$']}
            else:
                electrolyte = {'$z$': [self.z], r'$\epsilon_r$': [self.er0]}
            df_electrolyte = pd.DataFrame(data=electrolyte)
            df_pcsaft = pd.concat([df_pcsaft, df_electrolyte], axis=1)
        ref = {'reference': [self.reference]}
        df_ref = pd.DataFrame(data=ref)
        df_pcsaft = pd.concat([df_pcsaft, df_ref], axis=1)
        
        display(HTML(df_pcsaft.to_html(index=False)))
   

class mixture(object):
    '''
    class mixture
    Creates an object that cointains info about a mixture.

    Parameters
    ----------
    component1 : object
        component created with component class
    component2 : object
        component created with component class

    Attributes
    ----------
    name : list
        Name of the component
    ms : list
        list of chain lenght
    sigma : list
        list of segment diameter [m]
    eps : list
        List of dispersion energy energy [J]
    eAB : list
        List of Association Energy [J]
    kappaAB : list
        List of Association volumen
    sites = list
        Association sites [Bipolar, Positive, Negative]
    mupol : list
        List of dipolar moment [Debye]
    xpol : list
        List of fraction of dipolar segment on a chain molecule JC polar term
    z : float
        List of ions valence
    er : float
        List of dielectric constant
    cii : list
        polynomial coefficient for influence parameter in SGT [J m^5/mol^2]
    
    
    Methods
    -------
    add_component : adds a component to the mixture
    copy: returns a copy of the object
    kij_saft : add kij matrix for PC-SAFT
    kepsij_saft : add kepsij matrix for PC-SAFT
    lij_saft : add lij matrix for PC-SAFT
    ci : computes cij matrix at T for SGT (for SGTpy module)
    '''

    def __init__(self, component1, component2):
        self.components = [component1, component2]
        self.names = [component1.name, component2.name]
        self.cii = [component1.cii, component2.cii]
        self.Mw = [component1.Mw, component2.Mw]
        self.nc = 2

        self.sigmaT_bool = component1.sigmaT_bool + component2.sigmaT_bool > 0
        
 
        self.sigma0 = [component1.sigma0, component2.sigma0]
        self.t1 = [component1.t1, component2.t1]
        self.t2 = [component1.t2, component2.t2]
        self.t3 = [component1.t3, component2.t3]
        self.t4 = [component1.t4, component2.t4]
 
        self.eps = [component1.eps, component2.eps]
        self.ms = [component1.ms, component2.ms]
        

        self.eAB = [component1.eAB, component2.eAB]
        self.kappaAB = [component1.kappaAB, component2.kappaAB]
        self.sitesmix = [component1.sites, component2.sites]

        self.mupol = [component1.mupol, component2.mupol]
        self.xpol = [component1.xpol, component2.xpol]
        
        self.z =  [component1.z, component2.z]
        self.er0 = [component1.er0, component2.er0]
        self.er1 = [component1.er1, component2.er1]
        self.er2 = [component1.er2, component2.er2]

        self.viscosity_parameters = [component1.viscosity_parameters, component2.viscosity_parameters]
        
        
        ## kij matrix
        self.KIJ0saft = np.zeros([self.nc, self.nc])
        self.KIJ1saft = np.zeros([self.nc, self.nc])
        self.KIJ2saft = np.zeros([self.nc, self.nc])
        self.KIJ3saft = np.zeros([self.nc, self.nc])
        
        ## kij matrix
        self.KepsIJ0saft = np.zeros([self.nc, self.nc])
        self.KepsIJ1saft = np.zeros([self.nc, self.nc])
        self.KepsIJ2saft = np.zeros([self.nc, self.nc])
        self.KepsIJ3saft = np.zeros([self.nc, self.nc])

        ## lij matrix
        self.LIJ0saft = np.zeros([self.nc, self.nc])
        self.LIJ1saft = np.zeros([self.nc, self.nc])
        self.LIJ2saft = np.zeros([self.nc, self.nc])
        self.LIJ3saft = np.zeros([self.nc, self.nc])

        self.reference = [component1.reference, component2.reference]

    def add_component(self, component):
        """
        add_component(component)

        Method that add a component to the mixture

        Parameters
        ----------
        component : object
            pure fluid created with component function
        """
        self.components.append(component)
        self.names.append(component.name)
        self.cii.append(component.cii)
        self.Mw.append(component.Mw)

      
        self.sigmaT_bool = self.sigmaT_bool + component.sigmaT_bool > 0
                
         
        self.sigma0.append(component.sigma0)
        self.t1.append(component.t1)
        self.t2.append(component.t2)
        self.t3.append(component.t3)
        self.t4.append(component.t4)
                
        self.eps.append(component.eps)
        self.ms.append(component.ms)
        
        self.eAB.append(component.eAB)
        self.kappaAB.append(component.kappaAB)
        self.sitesmix.append(component.sites)

        self.mupol.append(component.mupol)
        self.xpol.append(component.xpol)

        self.z.append(component.z)
        self.er0.append(component.er0)
        self.er1.append(component.er1)
        self.er2.append(component.er2)

        self.viscosity_parameters.append(component.viscosity_parameters)

        self.nc += 1

        # kij matrix
        self.KIJ0saft = np.zeros([self.nc, self.nc])
        self.KIJ1saft = np.zeros([self.nc, self.nc])
        self.KIJ2saft = np.zeros([self.nc, self.nc])
        self.KIJ3saft = np.zeros([self.nc, self.nc])
        
        # kepsij matrix
        self.KepsIJ0saft = np.zeros([self.nc, self.nc])
        self.KepsIJ1saft = np.zeros([self.nc, self.nc])
        self.KepsIJ2saft = np.zeros([self.nc, self.nc])
        self.KepsIJ3saft = np.zeros([self.nc, self.nc])

        # lij matrix
        self.LIJ0saft = np.zeros([self.nc, self.nc])
        self.LIJ1saft = np.zeros([self.nc, self.nc])
        self.LIJ2saft = np.zeros([self.nc, self.nc])
        self.LIJ3saft = np.zeros([self.nc, self.nc])

        self.reference.append(component.reference)

    def __add__(self, new_component):
        if isinstance(new_component, component):
            self.add_component(new_component)
        else:
            raise Exception('You can only add components objects to an existing mixture')
        return self

    def kij_saft(self, kij0, kij1=None, kij2=None, kij3=None):
        r"""
        kij_saft(kij)

        Method that adds kij correction for energy dispersion 

        .. math::
            \epsilon_{ij} = (1 - k_{ij}) \sqrt{\epsilon_i \epsilon_j}

        kij correction is computed as follows:
        .. math::
            k_{ij} = k_{ij,0} + k_{ij,1} \cdot T +  k_{ij,2} \cdot T^2 + k_{ij,3} / T

        Parameters
        ----------
        kij0 : array_like
            kij0 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [Adim]
        kij1 : array_like, optional
            kij1 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K].
            If None, then a zero matrix is assumed.
        kij2 : array_like, optional
            kij2 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K^2].
            If None, then a zero matrix is assumed.
        kij3 : array_like, optional
            kij3 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [K].
            If None, then a zero matrix is assumed.

        """
        nc = self.nc
        KIJ0 = np.asarray(kij0)
        shape = KIJ0.shape

        isSquare = shape == (nc, nc)
        isSymmetric = np.allclose(KIJ0, KIJ0.T)
        diagZero = np.all(np.diagonal(KIJ0) == 0.)

        if isSquare and isSymmetric and diagZero:
            self.KIJ0saft = KIJ0
        else:
            raise Exception('kij0 matrix is not square, symmetric or diagonal==0')

        if kij1 is None:
            KIJ1 = np.zeros_like(kij0)
            self.KIJ1saft = KIJ1
        else:
            KIJ1 = np.asarray(kij1)
            shape = KIJ1.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(KIJ1, KIJ1.T)
            diagZero = np.all(np.diagonal(KIJ1) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.KIJ1saft = KIJ1
            else:
                raise Exception('kij1 matrix is not square, symmetric or diagonal==0')

        if kij2 is None:
            KIJ2 = np.zeros_like(kij0)
            self.KIJ2saft = KIJ2
        else:
            KIJ2 = np.asarray(kij2)
            shape = KIJ2.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(KIJ2, KIJ2.T)
            diagZero = np.all(np.diagonal(KIJ2) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.KIJ2saft = KIJ2
            else:
                raise Exception('kij2 matrix is not square, symmetric or diagonal==0')

        if kij3 is None:
            KIJ3 = np.zeros_like(kij0)
            self.KIJ3saft = KIJ3
        else:
            KIJ3 = np.asarray(kij3)
            shape = KIJ3.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(KIJ3, KIJ3.T)
            diagZero = np.all(np.diagonal(KIJ3) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.KIJ3saft = KIJ3
            else:
                raise Exception('kij3 matrix is not square, symmetric or diagonal==0')



    def kepsij_saft(self, kepsij0, kepsij1=None, kepsij2=None, kepsij3=None):
        r"""
        kepsij_saft(kij)

        Method that adds kij correction for energy dispersion 

        .. math::
            \epsilon_{AiBj} = (1 - keps_{ij}) \sqrt{\epsilon_{Ai} \epsilon_{Bj}}

        kepsij correction is computed as follows:
        .. math::
            keps_{ij} = keps_{ij,0} + keps_{ij,1} \cdot T +  keps_{ij,2} \cdot T^2 + keps_{ij,3} / T

        Parameters
        ----------
        kepsij0 : array_like
            kepsij0 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [Adim]
        kepsij1 : array_like, optional
            kepsij1 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K].
            If None, then a zero matrix is assumed.
        kepsij2 : array_like, optional
            kepsij2 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K^2].
            If None, then a zero matrix is assumed.
        kepsij3 : array_like, optional
            kepsij3 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [K].
            If None, then a zero matrix is assumed.

        """
        nc = self.nc
        KepsIJ0 = np.asarray(kepsij0)
        shape = KepsIJ0.shape

        isSquare = shape == (nc, nc)
        isSymmetric = np.allclose(KepsIJ0, KepsIJ0.T)
        diagZero = np.all(np.diagonal(KepsIJ0) == 0.)

        if isSquare and isSymmetric and diagZero:
            self.KepsIJ0saft = KepsIJ0
        else:
            raise Exception('kepsij0 matrix is not square, symmetric or diagonal==0')

        if kepsij1 is None:
            KepsIJ1 = np.zeros_like(kepsij0)
            self.KepsIJ1saft = KepsIJ1
        else:
            KepsIJ1 = np.asarray(kepsij1)
            shape = KepsIJ1.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(KepsIJ1, KepsIJ1.T)
            diagZero = np.all(np.diagonal(KepsIJ1) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.KepsIJ1saft = KepsIJ1
            else:
                raise Exception('kepsij1 matrix is not square, symmetric or diagonal==0')

        if kepsij2 is None:
            KepsIJ2 = np.zeros_like(kepsij0)
            self.KepsIJ2saft = KepsIJ2
        else:
            KepsIJ2 = np.asarray(kepsij2)
            shape = KepsIJ2.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(KepsIJ2, KepsIJ2.T)
            diagZero = np.all(np.diagonal(KepsIJ2) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.KepsIJ2saft = KepsIJ2
            else:
                raise Exception('kepsij2 matrix is not square, symmetric or diagonal==0')

        if kepsij3 is None:
            KepsIJ3 = np.zeros_like(kepsij0)
            self.KepsIJ3saft = KepsIJ3
        else:
            KepsIJ3 = np.asarray(kepsij3)
            shape = KepsIJ3.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(KepsIJ3, KepsIJ3.T)
            diagZero = np.all(np.diagonal(KepsIJ3) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.KepsIJ3saft = KepsIJ3
            else:
                raise Exception('kepsij3 matrix is not square, symmetric or diagonal==0')



    def lij_saft(self, lij0, lij1=None, lij2=None, lij3=None):
        r"""
        lij_saft(lij)

        Method that adds lij correction for segment diameter

        .. math::
            \sigma_{ij} = (1 - l_{ij})\frac{\sigma_i + \sigma_j}{2}

        lij correction is computed as follows:
        .. math::
            l_{ij} = l_{ij,0} + l_{ij,1} \cdot T +  l_{ij,2} \cdot T^2 + l_{ij,3} / T

        Parameters
        ----------
        lij0 : array_like
            lij0 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [Adim]
        lij1 : array_like, optional
            lij1 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K].
            If None, then a zero matrix is assumed.
        lij2 : array_like, optional
            lij2 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [1/K^2].
            If None, then a zero matrix is assumed.
        lij3 : array_like, optional
            lij3 matrix (Symmetric, Diagonal==0, shape=(nc, nc)) [K]
            If None, then a zero matrix is assumed.

        """

        nc = self.nc
        LIJ0 = np.asarray(lij0)
        shape = LIJ0.shape

        isSquare = shape == (nc, nc)
        isSymmetric = np.allclose(LIJ0, LIJ0.T)
        diagZero = np.all(np.diagonal(LIJ0) == 0.)

        if isSquare and isSymmetric and diagZero:
            self.LIJ0saft = LIJ0
        else:
            raise Exception('lij0 matrix is not square, symmetric or digonal==0')

        if lij1 is None:
            LIJ1 = np.zeros_like(lij0)
            self.LIJ1saft = LIJ1
        else:
            LIJ1 = np.asarray(lij1)
            shape = LIJ1.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(LIJ1, LIJ1.T)
            diagZero = np.all(np.diagonal(LIJ1) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.LIJ1saft = LIJ1
            else:
                raise Exception('lij1 matrix is not square, symmetric or diagonal==0')

        if lij2 is None:
            LIJ2 = np.zeros_like(lij0)
            self.LIJ2saft = LIJ2
        else:
            LIJ2 = np.asarray(lij2)
            shape = LIJ2.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(LIJ2, LIJ2.T)
            diagZero = np.all(np.diagonal(LIJ2) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.LIJ2saft = LIJ2
            else:
                raise Exception('lij2 matrix is not square, symmetric or diagonal==0')

        if lij3 is None:
            LIJ3 = np.zeros_like(lij0)
            self.LIJ3saft = LIJ3
        else:
            LIJ3 = np.asarray(lij3)
            shape = LIJ3.shape

            isSquare = shape == (nc, nc)
            isSymmetric = np.allclose(LIJ3, LIJ3.T)
            diagZero = np.all(np.diagonal(LIJ3) == 0.)

            if isSquare and isSymmetric and diagZero:
                self.LIJ3saft = LIJ3
            else:
                raise Exception('lij3 matrix is not square, symmetric or diagonal==0')

    def set_kijsaft(self, i, j, kij0, kij1=0., kij2=0., kij3=0.):
        r"""
        set_kijsaft(i,j, kij0, kij1, kij2, kij3)

        Method that sets the kij correction for dispersion energy
        between component i and component j.

        .. math::
            \epsilon_{AiBj} = (1 - keps_{ij}) \dfrac{\epsilon_{Ai} + \epsilon_{Bj}}{2}
            
        kij correction is computed as follows:
        .. math::
            k_{ij} = k_{ij,0} + k_{ij,1} \cdot T +  k_{ij,2} \cdot T^2 + k_{ij,3} / T

        Parameters
        ----------
        i : int
            index of component i.
        j : int
            index of component j.
        kij0 : float
            kij0 value between component i and j [Adim]
        kij1 : float, optional
            kij1 value between component i and j [1/K]. Default to zero.
        kij2 : float, optional
            kij2 value between component i and j [1/K^2]. Default to zero.
        kij3 : float, optional
            kij3 value between component i and j [K]. Default to zero.

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
            raise Exception('Cannot set kij for i=j')

        if typei and typej and nc_i and nc_j and i_j:
            self.KIJ0saft[i, j] = kij0
            self.KIJ0saft[j, i] = kij0

            self.KIJ1saft[i, j] = kij1
            self.KIJ1saft[j, i] = kij1

            self.KIJ2saft[i, j] = kij2
            self.KIJ2saft[j, i] = kij2

            self.KIJ3saft[i, j] = kij3
            self.KIJ3saft[j, i] = kij3

            
            
    def set_kepsijsaft(self, i, j, kepsij0, kepsij1=0., kepsij2=0., kepsij3=0.):
         r"""
         set_kepsijsaft(i, j, kepsij0, kepsij1, kepsij2, kepsij3)

         Method that sets the kij correction for association energy
         between component i and component j.

         .. math::
             \epsilon_{AiBj} = (1 - keps_{ij}) \sqrt{\epsilon_{Ai} \epsilon_{Bj}}

         kij correction is computed as follows:
         .. math::
             keps_{ij} = keps_{ij,0} + keps_{ij,1} \cdot T +  keps_{ij,2} \cdot T^2 + keps_{ij,3} / T

         Parameters
         ----------
         i : int
             index of component i.
         j : int
             index of component j.
         kepsij0 : float
             kepsij0 value between component i and j [Adim]
         kepsij1 : float, optional
             kepsij1 value between component i and j [1/K]. Default to zero.
         kepsij2 : float, optional
             kepsij2 value between component i and j [1/K^2]. Default to zero.
         kepsij3 : float, optional
             kepsij3 value between component i and j [K]. Default to zero.

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
             raise Exception('Cannot set kepsij for i=j')

         if typei and typej and nc_i and nc_j and i_j:
             self.KepsIJ0saft[i, j] = kepsij0
             self.KepsIJ0saft[j, i] = kepsij0

             self.KepsIJ1saft[i, j] = kepsij1
             self.KepsIJ1saft[j, i] = kepsij1

             self.KepsIJ2saft[i, j] = kepsij2
             self.KepsIJ2saft[j, i] = kepsij2

             self.KepsIJ3saft[i, j] = kepsij3
             self.KepsIJ3saft[j, i] = kepsij3

    def set_lijsaft(self, i, j, lij0, lij1=0., lij2=0., lij3=0.):
        r"""
        set_lijsaft(i,j, lij0, lij1, lij2, lij3)

        Method that sets the lij correction for segment diameter
        between component i and component j.

        .. math::
            \sigma_{ij} = (1 - l_{ij})\frac{\sigma_i + \sigma_j}{2}

        lij correction is computed as follows:
        .. math::
            l_{ij} = l_{ij,0} + l_{ij,1} \cdot T +  l_{ij,2} \cdot T^2 + l_{ij,3} / T

        Parameters
        ----------
        i : int
            index of component i.
        j : int
            index of component j.
        lij0 : float
            lij0 value between component i and j [Adim]
        lij1 : float, optional
            lij1 value between component i and j [1/K]. Default to zero.
        lij2 : float, optional
            lij2 value between component i and j [1/K^2]. Default to zero.
        lij3 : float, optional
            lij3 value between component i and j [K]. Default to zero.

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
            raise Exception('Cannot set lij for i=j')

        if typei and typej and nc_i and nc_j and i_j:
            self.LIJ0saft[i, j] = lij0
            self.LIJ0saft[j, i] = lij0

            self.LIJ1saft[i, j] = lij1
            self.LIJ1saft[j, i] = lij1

            self.LIJ2saft[i, j] = lij2
            self.LIJ2saft[j, i] = lij2

            self.LIJ3saft[i, j] = lij3
            self.LIJ3saft[j, i] = lij3

    def ci(self, T):
        """
        Method that computes the matrix of cij interaction parameter for SGT at
        given temperature (to use with sgtpy module).

        Parameters
        ----------
        T : float
            absolute temperature [K]

        Returns
        -------
        ci : array_like
            influence parameter matrix at given temperature [J m^5 / mol^2]

        """

        n = len(self.cii)
        ci = np.zeros(n)
        for i in range(n):
            ci[i] = np.polyval(self.cii[i], T)
        self.cij = np.sqrt(np.outer(ci, ci))
        return self.cij

    def copy(self):
        """
        Method that return a copy of the mixture

        Returns
        -------
        mix : object
            returns a copy a of the mixture
        """
        return copy(self)
    
    def printParameters(self):
        """
        Method that prints the PC-SAFT parameters for a mixture.
        """
        sigmaT = []
        erT = []
        for i in range(self.nc):
            erT.append(str(self.er0[i]) + ' + ' + str(self.er1[i]) + 'T + ' + str(self.er2[i]) + 'ln T')
            if self.t1[i] != 0:
                sigmaT.append('$\sigma(T)$')
            else:
                sigmaT.append(self.sigma0[i]*1e10 )
        pcsaft = {'component': self.names, 'Mw': self.Mw, '$ms$': self.ms, 
                  '$\sigma$': sigmaT, r'$\epsilon / k_B $': [i/kb for i in self.eps]}
        df_pcsaft = pd.DataFrame(data=pcsaft)
        if max(max(self.sitesmix)) != 0:
            assoc = {'[B, P, N]': self.sitesmix, r'$\epsilon _{AB} / k_B$': [i/kb for i in self.eAB], r'$k_{AB}$': self.kappaAB}
            df_assoc = pd.DataFrame(data=assoc)
            df_pcsaft = pd.concat([df_pcsaft, df_assoc], axis=1)
        if max(self.xpol) != 0:
            polar = {'$x_p$': self.xpol, r'$\mu$': self.mupol}
            polar = pd.DataFrame(data=polar)
            df_pcsaft = pd.concat([df_pcsaft, polar], axis=1)
        if max(self.z) != 0:
            electrolyte = {'$z$': self.z, r'$\epsilon_r(T)$': erT}
            electrolyte = pd.DataFrame(data=electrolyte)
            df_pcsaft = pd.concat([df_pcsaft, electrolyte], axis=1)
        ref = {'reference': self.reference}
        df_ref = pd.DataFrame(data=ref)
        df_pcsaft = pd.concat([df_pcsaft, df_ref], axis=1)
        display(HTML(df_pcsaft.to_html(index=False)))

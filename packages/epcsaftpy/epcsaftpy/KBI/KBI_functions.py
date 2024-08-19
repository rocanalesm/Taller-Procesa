import numpy as np


kB = 1.3806488e-23 # J K-1
Na = 6.022142e23   # mol -1
R = kB*Na

# Function to save information as a dictonary

class Result(dict):

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)


    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):

        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in self.items()])
        else:
            return self.__class__.__name__ + "()"
        
# Function to compute numerical derivates
def derivate(fun, x, temp_aux, P, saft, rho0, Xass0):
    nc = saft.nc
    f = fun(x, temp_aux, P, saft, rho0, Xass0)
    n = np.size(f)
    h = np.finfo(1.0).eps
    if n>1:
        df = np.zeros([nc, n])
    else:
        df = np.zeros(nc)
        
    for i in range(nc):
        dx = np.zeros(nc, dtype = 'complex128')
        try:
            Xass0n = np.zeros_like(Xass0, dtype = 'complex128')
            Xass0n[:] = Xass0
        except:
            Xass0n = None

        dx[i] = h * 1j     
        out = fun(x + dx, temp_aux, P, saft, rho0, Xass0n)
        df[i] = out.imag/h
    return df.T        

       
# Function to compute volume at T and P constant.
def V(N, temp_aux, P, saft, rho0, Xass0):
    NT = np.sum(N)
    x = N/NT
    rho, Xass0 = saft.density_aux(x, temp_aux, P, "L", rho0, Xass0)
    return NT/rho
    
# Function to compute chemical potential at T and P constant.
def mu(N, temp_aux, P, saft, rho0, Xass0):
    NT = np.sum(N)
    x = N/NT
    rho, Xass0 = saft.density_aux(x, temp_aux, P, "L", rho0, Xass0)
    mu, _ = saft.muad_aux(rho*x, temp_aux, Xass0)
    return mu

# Function to compute the isothermal compressibility
def kTfun(x, temp_aux, rho, saft, Xass0):  
    _, dPsaft, _ = saft.dP_drho_aux(x, rho, temp_aux, Xass0)
    return 1/dPsaft/rho 
    
    

# Function to compute KBI�s

def KBI_binary(x, T, P, saft, rho0 = None, Xass0 = None , full_output = False):
    

    temp_aux  = saft.temperature_aux(T)
    
    rho, Xass0 = saft.density_aux(x, temp_aux, P, "L", rho0, Xass0)
    rhoi = x*rho
    
    # Partial molar volume (array)
    vmpi = derivate(V, x, temp_aux, P, saft, rho0, Xass0) 
    
    # muij
    muij = derivate(mu, x, temp_aux, P, saft, rho0, Xass0) * x
    
    # kT 
    kT = kTfun(x, temp_aux, rho, saft, Xass0)
    
    phi = rhoi*vmpi
    
    


    N11 = -1 + rhoi[0]*R*T*kT + phi[1]**2/muij[0,0]
    N22 = -1 + rhoi[1]*R*T*kT + phi[0]**2/muij[1,1]
    N12 = rhoi[1]*R*T*kT - phi[0]*phi[1]/muij[1,1]

    G11 = N11/rhoi[0]
    G22 = N22/rhoi[1]
    G12 = N12/rhoi[1]
    
    Gij = np.array([[G11, G12],
                   [G12, G22]])
    

        
    if full_output:
        sol = {'T' : T, 'P': P, 'rho':rho, 'Xass':Xass0 , 'vmpi':vmpi,
        'muij' : muij, 'kT':kT, 'Gij' : Gij}
        out = Result(sol)
        return out    
    return Gij



# Function to compute KBI

def KBI(x, T, P, saft, rho0 = None, Xass0 = None , full_output = False):
    

    temp_aux  = saft.temperature_aux(T)
    nc = saft.nc
    
    rho, Xass0 = saft.density_aux(x, temp_aux, P, "L", rho0, Xass0)
    rhoi = x*rho
    
    # Partial molar volume (array)
    vmpi = derivate(V, x, temp_aux, P, saft, rho0, Xass0) 
    
    # muij
    muij = derivate(mu, x, temp_aux, P, saft, rho0, Xass0) * x
    muij  = muij.T
    
    # kT 
    kT = kTfun(x, temp_aux, rho, saft, Xass0)
    
    phii = rhoi*vmpi
    muij_vmpi = np.append(muij, [vmpi], axis=0)
    phi_kt = np.append(-phii, [R*T*kT], axis=0) 

    k = 0
    Gij = np.zeros([nc, nc])
    for k in range(nc):
        Mki = np.delete(muij_vmpi, k, 0)
        bk = np.delete(phi_kt, k, 0) 
        dki = np.zeros(nc)
        dki[k] = 1
        Gij[k] = (np.linalg.solve(Mki, bk) - dki)/rhoi
    Gij
    

        
    if full_output:
        sol = {'T' : T, 'P': P, 'rho':rho, 'Xass':Xass0 , 'vmpi':vmpi,
        'muij' : muij, 'kT':kT, 'Gij' : Gij}
        out = Result(sol)
        return out    
    return Gij

def TH(x, T, P, saft, rho0 = None, Xass0 = None , full_output = False):
    temp_aux  = saft.temperature_aux(T)
    rho, Xass0 = saft.density_aux(x, temp_aux, P, "L", rho0, Xass0)
    # chemical potential derivative (matrix)
    muij = derivate(mu, x, temp_aux, P, saft, rho0, Xass0) * x
    muij  = muij.T
    nc=saft.nc
    nTH=saft.nc-1
    THij = muij[:nTH,:nTH]
    massbalancecorrection=np.vstack([muij[:nTH,-1]]*nTH)
    THij -= massbalancecorrection
    if full_output:
        sol = {'T' : T, 'P': P, 'rho':rho, 'Xass':Xass0 , 
        'muij' : muij, 'THij' : THij}
        out = Result(sol)
        return out    
    return np.squeeze(THij)

def THdiag(x, T, P, saft, rho0 = None, Xass0 = None , full_output = False):
    temp_aux  = saft.temperature_aux(T)
    rho, Xass0 = saft.density_aux(x, temp_aux, P, "L", rho0, Xass0)
    # chemical potential derivative (matrix)
    muij = derivate(mu, x, temp_aux, P, saft, rho0, Xass0) * x
    muij  = muij.T
    nc=saft.nc
    nTH=saft.nc-1
    THij1 = muij[:nTH,:nTH]
    massbalancecorrection1=np.vstack([muij[:nTH,-1]]*nTH)
    THij1 -= massbalancecorrection1
    THij2 = muij[1:,1:]
    massbalancecorrection2=np.vstack([muij[1:,0]]*nTH)
    THij2 -= massbalancecorrection2
    THii= np.hstack((np.diag(THij1),np.diag(THij2)[1:]))

    if full_output:
        sol = {'T' : T, 'P': P, 'rho':rho, 'Xass':Xass0 , 
        'muij' : muij, 'THij' : THdiag}
        out = Result(sol)
        return out    
    return np.squeeze(THii)
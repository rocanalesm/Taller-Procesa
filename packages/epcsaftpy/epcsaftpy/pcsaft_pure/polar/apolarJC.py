from __future__ import division, print_function, absolute_import
import numpy as np

def apolarJC(eta, ms, dia3, betaxpmu2):
    I2p = (625*(58212*eta**3 - 28845*eta**2*np.pi - 5427*eta*np.pi**2 + 2500*np.pi**3))/((3927*eta - 1250*np.pi)**2*np.pi)
    
    I3p = 1 + (9*eta*(-31717*eta + 20239*np.pi))/(180531*eta**2 - 88584*eta*np.pi + 25000*np.pi**2)
    
    a2 = (-4*betaxpmu2**2*eta*I2p*ms)/(3.*dia3**2)
    
    a3 = (10*betaxpmu2**3*eta**2*I3p*ms)/(9.*dia3**3)

    a = a2/(1 - a3/a2)
    return a

def dapolarJC_deta(eta, ms, dia3, betaxpmu2):
    I2p = (625*(58212*eta**3 - 28845*eta**2*np.pi - 5427*eta*np.pi**2 + 2500*np.pi**3))/((3927*eta - 1250*np.pi)**2*np.pi)
    dI2p = (1875*(76199508*eta**3 - 72765000*eta**2*np.pi + 31141443*eta*np.pi**2 - 4283750*np.pi**3))/((3927*eta - 1250*np.pi)**3*np.pi)
    
    I3p = 1 + (9*eta*(-31717*eta + 20239*np.pi))/(180531*eta**2 - 88584*eta*np.pi + 25000*np.pi**2)
    dI3p = (9*np.pi*(-844148181*eta**2 - 1585850000*eta*np.pi + 505975000*np.pi**2))/(180531*eta**2 - 88584*eta*np.pi + 25000*np.pi**2)**2
    
    a2 = (-4*betaxpmu2**2*eta*I2p*ms)/(3.*dia3**2)
    da2 = a2*(1/eta + dI2p/I2p)
    
    a3 = (10*betaxpmu2**3*eta**2*I3p*ms)/(9.*dia3**3)
    da3 = a3*(2/eta + dI3p/I3p)

    a = a2/(1 - a3/a2)
    da = (a2*((a2 - 2*a3)*da2 + a2*da3))/(a2 - a3)**2
    return np.array([a, da])


def d2apolarJC_deta(eta, ms, dia3, betaxpmu2):
    I2p = (625*(58212*eta**3 - 28845*eta**2*np.pi - 5427*eta*np.pi**2 + 2500*np.pi**3))/((3927*eta - 1250*np.pi)**2*np.pi)
    dI2p = (1875*(76199508*eta**3 - 72765000*eta**2*np.pi + 31141443*eta*np.pi**2 - 4283750*np.pi**3))/((3927*eta - 1250*np.pi)**3*np.pi)
    d2I2p = (-11250*(10445398887*eta - 1923342500*np.pi)*np.pi)/(3927*eta - 1250*np.pi)**4
    
    I3p = 1 + (9*eta*(-31717*eta + 20239*np.pi))/(180531*eta**2 - 88584*eta*np.pi + 25000*np.pi**2)
    dI3p = (9*np.pi*(-844148181*eta**2 - 1585850000*eta*np.pi + 505975000*np.pi**2))/(180531*eta**2 - 88584*eta*np.pi + 25000*np.pi**2)**2
    d2I3p = (18*np.pi*(152394915264111*eta**3 + 429442629525000*eta**2*np.pi - 274032518175000*eta*np.pi**2 + 24998164400000*np.pi**3))/(180531*eta**2 - 88584*eta*np.pi + 25000*np.pi**2)**3
    
    a2 = (-4*betaxpmu2**2*eta*I2p*ms)/(3.*dia3**2)
    da2 = a2*(1/eta + dI2p/I2p)
    d2a2 = a2*((2*dI2p + d2I2p*eta)/(eta*I2p))
    
    a3 = (10*betaxpmu2**3*eta**2*I3p*ms)/(9.*dia3**3)
    da3 = a3*(2/eta + dI3p/I3p)
    d2a3 = a3*((4*dI3p*eta + d2I3p*eta**2 + 2*I3p)/(eta**2*I3p))

    a = a2/(1 - a3/a2)
    da = (a2*((a2 - 2*a3)*da2 + a2*da3))/(a2 - a3)**2
    d2a = (a2**3*(d2a2 + d2a3) + 2*a3**2*da2**2 + 2*a2*a3*(a3*d2a2 - 2*da2*da3) - a2**2*(a3*(3*d2a2 + d2a3) - 2*da3**2))/(a2 - a3)**3
    return np.array([a, da, d2a])

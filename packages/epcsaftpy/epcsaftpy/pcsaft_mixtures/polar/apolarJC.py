from __future__ import division, print_function, absolute_import
import numpy as np


def apolarJC(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3):
    xij = np.multiply.outer(x, x)
    xijk = np.multiply.outer(x, xij)
    
    auxa2 = np.sum(xij*aux1_polarij/dij3) / 9
    if auxa2==0:
        a = 0.
    else:    
        rhoad = eta/(np.pi/6)
        rhoad2, rhoad3 = rhoad**2, rhoad**3
        I2p = 1 - 0.3618*rhoad - 0.3205*rhoad2 + 0.1078*rhoad3
        I2p /= (1 - 0.5236*rhoad)**2
        I3p = 1 + 0.62378*rhoad - 0.11658*rhoad2
        I3p /= 1 - 0.59056*rhoad + 0.20059*rhoad2
        
        pirho = np.pi * rho
        a2 = -2 * I2p * pirho
        a2 *= auxa2
        
        a3 = 5 * I3p * pirho**2
        a3 *= np.sum(xijk*aux2_polarijk/dijk3) / 162  
    
        a = a2/(1 - a3/a2)
    
    return a


def dapolarJC_drho(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3):
    xij = np.multiply.outer(x, x)
    xijk = np.multiply.outer(x, xij)
    
    auxa2 = -2 * np.sum(xij*aux1_polarij/dij3) / 9
    if auxa2 == 0:
        a, da = 0., 0.
    else:  
        rhoad = eta/(np.pi/6)
        rhoad2, rhoad3 = rhoad**2, rhoad**3
        drhoad_drho = rhoad/rho
        
        aux2p = 1 - 0.5236*rhoad
        I2p = 1 - 0.3618*rhoad - 0.3205*rhoad2 + 0.1078*rhoad3
        I2p /= aux2p**2
        dI2p = 0.6854 - 0.83043848*rhoad 
        dI2p += 0.3234*rhoad2 - 0.05644408*rhoad3
        dI2p /= aux2p**3
        dI2p *= drhoad_drho
        
        aux3p = 1 - 0.59056*rhoad + 0.20059*rhoad2
        I3p = 1 + 0.62378*rhoad - 0.11658*rhoad2
        I3p /= aux3p
        dI3p = 1.21434 - 0.63434*rhoad
        dI3p -= 0.0562765454*rhoad2
        dI3p /= aux3p**2
        dI3p *= drhoad_drho
        
        pirho = np.pi * rho
        a2 = I2p * pirho
        a2 *= auxa2
        da2 = I2p * np.pi 
        da2 +=  dI2p * pirho
        da2 *= auxa2
        
    
        pirho2 = pirho * pirho
        auxa3 = 5 * np.sum(xijk*aux2_polarijk/dijk3) / 162  
        a3 = I3p * pirho2
        a3 *= auxa3 
        da3 = 2 * I3p * pirho2 / rho
        da3 += dI3p * pirho2
        da3 *= auxa3 
    
        a3a2 = a3/a2
        auxa = 1 - a3a2
        a = a2/auxa
        da = (auxa - a3a2)*da2 + da3
        da /= auxa**2
    
    return np.array([a, da])

def d2apolarJC_drho(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3):
    xij = np.multiply.outer(x, x)
    xijk = np.multiply.outer(x, xij)
    auxa2 = np.sum(xij*aux1_polarij/dij3) / 9
    if auxa2 == 0:
        a, da, d2a = 0., 0., 0.
    else:
        rhoad = eta/(np.pi/6)
        rhoad2, rhoad3 = rhoad**2, rhoad**3
        drhoad_drho = rhoad/rho
        drhoad_drho2 = drhoad_drho**2
        
        aux2p = 1 - 0.5236*rhoad
        I2p = 1 - 0.3618*rhoad - 0.3205*rhoad2 + 0.1078*rhoad3
        I2p /= aux2p**2
        dI2p = 0.6854 - 0.83043848*rhoad 
        dI2p += 0.3234*rhoad2 - 0.05644408*rhoad3
        dI2p /= aux2p**3
        dI2p *= drhoad_drho
        d2I2p = 0.24618784000000002 - 0.222835176256*rhoad
        d2I2p /= aux2p**4
        d2I2p *= drhoad_drho2
        
        aux3p = 1 - 0.59056*rhoad + 0.20059*rhoad2
        I3p = 1 + 0.62378*rhoad - 0.11658*rhoad2
        I3p /= aux3p
        dI3p = 1.21434 - 0.63434*rhoad
        dI3p -= 0.0562765454*rhoad2
        dI3p /= aux3p**2
        dI3p *= drhoad_drho
        d2I3p = 0.7999412608 - 1.4615067636*rhoad
        d2I3p += 0.3817267818*rhoad2 + 0.022577024483572003*rhoad3
        d2I3p /= aux3p**3
        d2I3p *= drhoad_drho2
            
        pirho = np.pi * rho
        a2 = -2 * I2p * pirho
        a2 *= auxa2
        da2 = -2 * I2p * np.pi 
        da2 -=  2 * dI2p * pirho
        da2 *= auxa2
        d2a2 = -2 * dI2p * np.pi 
        d2a2 -=  2 * d2I2p * pirho
        d2a2 -=  2 * dI2p * np.pi
        d2a2 *= auxa2
        
    
        pirho2 = pirho * pirho
        auxa3 = np.sum(xijk*aux2_polarijk/dijk3) / 162  
        a3 = 5 * I3p * pirho2
        a3 *= auxa3 
        da3 = 10 * I3p * pirho2 / rho
        da3 += 5 * dI3p * pirho2
        da3 *= auxa3 
        d2a3 = dI3p * pirho2 / rho
        d2a3 += I3p * np.pi**2
        d2a3 *= 2
        d2a3 += d2I3p * pirho2
        d2a3 += 2 * dI3p * pirho2 / rho
        d2a3 *= 5 * auxa3 
    
        a3a2 = a3/a2
        aux1 = -da2/a2
        aux2 = da3/a3
        daux1 = d2a2/da2 - da2/a2
        daux1 *= aux1
        daux2 = d2a3/da3 - da3/a3
        daux2 *= aux2
        da3a2 = aux1 + aux2
        da3a2 *= a3a2
        d2a3a2 = (aux1 + aux2)*da3a2
        d2a3a2 += (daux1 + daux2)*a3a2
        auxa = 1 - a3a2
        auxa2 = auxa**2
        auxa3 = auxa2*auxa
        
        a = a2/auxa
        da = auxa*da2 + a2*da3a2
        da /= auxa2
        d2a = auxa2*d2a2 + 2*auxa*da2*da3a2
        d2a += a2*(auxa*d2a3a2 + 2*da3a2**2)
        d2a /= auxa3
    
    return np.array([a, da, d2a])



def d2apolarJC_dx(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3, deta_dx):
    xij = np.multiply.outer(x, x)
    xijk = np.multiply.outer(x, xij)
    auxij = aux1_polarij/dij3
    auxa2 = np.sum(xij*auxij)
    if auxa2 == 0:
        a, dax = 0., 0. * x
    else:    
        auxrhoad = np.pi/6
        rhoad = eta/auxrhoad
        rhoad2, rhoad3 = rhoad**2, rhoad**3
        drhoad_dx = deta_dx/auxrhoad
        
        aux2p = 1 - 0.5236*rhoad
        I2p = 1 - 0.3618*rhoad - 0.3205*rhoad2 + 0.1078*rhoad3
        I2p /= aux2p**2
        dI2px = 0.6854 - 0.83043848*rhoad 
        dI2px += 0.3234*rhoad2 - 0.05644408*rhoad3
        dI2px /= aux2p**3
        dI2px *= drhoad_dx
    
        
        aux3p = 1 - 0.59056*rhoad + 0.20059*rhoad2
        I3p = 1 + 0.62378*rhoad - 0.11658*rhoad2
        I3p /= aux3p
        dI3px = 1.21434 - 0.63434*rhoad
        dI3px -= 0.0562765454*rhoad2
        dI3px /= aux3p**2
        dI3px *= drhoad_dx
    
        pirho = np.pi * rho
        aux_a2 = -2 * pirho / 9
        a2 = I2p * auxa2
        a2 *= aux_a2
        da2x = dI2px * auxa2
        da2x += 2 * I2p * auxij@x
        da2x *= aux_a2
    
        
        auxijk = aux2_polarijk/dijk3
        auxa3 = np.sum(xijk *  auxijk)
        pirho2 = pirho * pirho
        aux_a3 = 5 * pirho2 / 162
        a3 = I3p * auxa3
        a3 *= aux_a3
        da3x = dI3px * auxa3
        da3x += 3 * I3p * np.sum(auxijk * xij, axis = (1, 2))
        da3x *= aux_a3 
     
    
        a3a2 = a3/a2
        auxa = (1 - a3a2)
        a = a2/auxa
        dax = (1 - 2*a3a2)*da2x + da3x
        dax /= auxa**2
           
    return a, dax


def d2apolarJC_dxrho(x, rho, eta, aux1_polarij, aux2_polarijk, dij3, dijk3, deta_dx):
    xij = np.multiply.outer(x, x)
    xijk = np.multiply.outer(x, xij)
    auxij = aux1_polarij/dij3
    auxa2 = np.sum(xij*auxij)
    if auxa2 == 0:
        a, da, dax = 0., 0., 0. * x
    else:
        auxrhoad = np.pi/6
        rhoad = eta/auxrhoad
        rhoad2, rhoad3 = rhoad**2, rhoad**3
        drhoad_dx = deta_dx/auxrhoad
        drhoad_drho = rhoad/rho
        
        aux2p = 1 - 0.5236*rhoad
        I2p = 1 - 0.3618*rhoad - 0.3205*rhoad2 + 0.1078*rhoad3
        I2p /= aux2p**2
        dI2p = 0.6854 - 0.83043848*rhoad 
        dI2p += 0.3234*rhoad2 - 0.05644408*rhoad3
        dI2p /= aux2p**3
        dI2p *= drhoad_drho
        dI2px = 0.6854 - 0.83043848*rhoad 
        dI2px += 0.3234*rhoad2 - 0.05644408*rhoad3
        dI2px /= aux2p**3
        dI2px *= drhoad_dx
    
        
        aux3p = 1 - 0.59056*rhoad + 0.20059*rhoad2
        I3p = 1 + 0.62378*rhoad - 0.11658*rhoad2
        I3p /= aux3p
        dI3p = 1.21434 - 0.63434*rhoad
        dI3p -= 0.0562765454*rhoad2
        dI3p /= aux3p**2
        dI3p *= drhoad_drho
        dI3px = 1.21434 - 0.63434*rhoad
        dI3px -= 0.0562765454*rhoad2
        dI3px /= aux3p**2
        dI3px *= drhoad_dx
    
        pirho = np.pi * rho
        aux_a2 = -2 * pirho / 9
        a2 = I2p * auxa2
        a2 *= aux_a2
        da2 = I2p / rho 
        da2 +=  dI2p
        da2 *= aux_a2
        da2x = dI2px * auxa2
        da2x += 2 * I2p * auxij@x
        da2x *= aux_a2
    
        auxijk = aux2_polarijk/dijk3
        auxa3 = np.sum(xijk *  auxijk)
        pirho2 = pirho * pirho
        aux_a3 = 5 * pirho2 / 162
        a3 = I3p * auxa3
        a3 *= aux_a3
        da3 = 2 * I3p / rho
        da3 += dI3p
        da3 *= aux_a3 
        da3x = dI3px * auxa3
        da3x += 3 * I3p * np.sum(auxijk * xij, axis = (1, 2))
        da3x *= aux_a3 
    
        a3a2 = a3/a2
        auxa = 1 - a3a2
        auxa2 = auxa**2
        auxa3a2 = auxa - a3a2
        a = a2/auxa
        da = auxa3a2 * da2 + da3
        da /= auxa2
        dax = auxa3a2 * da2x + da3x
        dax /= auxa2
       
    return np.array([a, da]), dax
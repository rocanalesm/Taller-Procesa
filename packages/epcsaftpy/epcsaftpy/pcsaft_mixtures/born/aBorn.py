from __future__ import division, print_function, absolute_import
import numpy as np

def aborn(xz2, er, dii, betae2_e0): 
    a = - np.sum(xz2/dii)
    a *= 1. - 1./er
    a *= betae2_e0/(4. * np.pi)
    return a

def daborn_drho(xz2, er, dii, betae2_e0): 
    a = - np.sum(xz2/dii)
    a *= 1. - 1./er
    a *= betae2_e0/(4. * np.pi)
    return np.array([a, 0.])

def d2aborn_drho(xz2, er, dii, betae2_e0): 
    a = - np.sum(xz2/dii)
    a *= 1. - 1./er
    a *= betae2_e0/(4. * np.pi)
    return np.array([a, 0., 0.])

def daborn_dx(xz2, z, er, derx, dii, betae2_e0): 
    aux_a = -betae2_e0/(4. * np.pi)
    sum_xz2 = np.sum(xz2/dii)
    dsum_xz2 = z**2/dii
    
    a = 1. - 1./er
    a *= np.sum(xz2/dii)
    a *= aux_a
    
    dax = sum_xz2 * derx
    dax += (er - 1) * er * dsum_xz2
    dax *= aux_a / er**2
    
    return a, dax

def daborn_dxrho(xz2, z, er, derx, dii, betae2_e0): 
    aux_a = -betae2_e0/(4. * np.pi)
    sum_xz2 = np.sum(xz2/dii)
    dsum_xz2 = z**2/dii
    
    a = 1. - 1./er
    a *= np.sum(xz2/dii)
    a *= aux_a
    
    dax = sum_xz2 * derx
    dax += (er - 1) * er * dsum_xz2
    dax *= aux_a / er**2
    
    return np.array([a, 0.]), dax


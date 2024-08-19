from __future__ import division, print_function, absolute_import
import numpy as np


def achain(x, ms, gHS, diag_index):  
    return - np.dot(x * (ms - 1), np.log(gHS[diag_index]))

def dachain_dxhi00(x, ms, gHS, diag_index):
    ghs, dghs = gHS
    ghsii = ghs[diag_index]
    dghsii = dghs[diag_index]
    
    lng = np.log(ghsii)
    
    dlng = dghsii
    dlng /= ghsii
    
    lngHS = np.array([lng, dlng])
    return - lngHS@(x * (ms - 1.))

def d2achain_dxhi00(x, ms, gHS, diag_index):
    ghs, dghs, d2ghs = gHS
    ghsii = ghs[diag_index]
    dghsii = dghs[diag_index]
    d2ghsii = d2ghs[diag_index]
    
    lng = np.log(ghsii)
    
    dlng = dghsii
    dlng /= ghsii
    
    d2lng = - dlng**2
    d2lng += d2ghsii/ghsii
    
    lngHS = np.array([lng, dlng, d2lng])
    return - lngHS@(x * (ms - 1.))

def dachain_dx(x, ms, gHS, dgHSx, diag_index):
    ghsii = gHS[diag_index]
    dgHSxii = np.einsum('ijj->ij',dgHSx) 
    
    lng = np.log(ghsii)
   
    
    dlngx = np.copy(dgHSxii)
    dlngx /= ghsii

    ms_1 = (ms - 1)
    aux_chain = x * ms_1
    a = - lng@aux_chain
    dax = -  dlngx@aux_chain
    dax -= ms_1*lng
    
    return a, dax


def dachain_dxxhi(x, ms, gHS, dgHSx, diag_index):
    ghs, dghs = gHS
    ghsii = ghs[diag_index]
    dghsii = dghs[diag_index]
    dgHSxii = np.einsum('ijj->ij',dgHSx) 
    
    lng = np.log(ghsii)
    
    dlng = dghsii
    dlng /= ghsii
    
    lngHS = np.array([lng, dlng])
    
    ms_1 = (ms - 1)
    aux_chain = x * ms_1
    a = - lngHS@aux_chain
    
    dlngx = np.copy(dgHSxii)
    dlngx /= ghsii
    dax = -  dlngx@aux_chain
    dax -= ms_1*lng
    
    return a, dax
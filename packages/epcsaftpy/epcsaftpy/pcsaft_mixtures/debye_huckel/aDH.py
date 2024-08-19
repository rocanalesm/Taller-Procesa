from __future__ import division, print_function, absolute_import
import numpy as np

def adh(rho, xz2, er, dii, betae2_e0): 
    betae2_e = betae2_e0 / er
    kappa = np.sqrt(betae2_e * rho * np.sum(xz2))
    
    kappaai = kappa *  dii
    kappaai_1 = 1 + kappaai
    
    aux_Xi1 = 1.5 - 2*kappaai_1 + kappaai_1**2/2. 
    aux_Xi1 += np.log(kappaai_1)
    with np.errstate(divide='ignore', invalid='ignore'):
        aux_Xi2 = 3/kappaai**3
        
        Xi = aux_Xi1 * aux_Xi2 
            
        sum_X = np.dot(Xi, xz2)
        
        aux_a = - betae2_e / (12 * np.pi)
        a = kappa * sum_X
        a *= aux_a
    return np.nan_to_num(a)


def dadh_drho(rho, xz2, er, dii, betae2_e0): 
    betae2_e = betae2_e0 / er
    kappa = np.sqrt(betae2_e * rho * np.sum(xz2))
    dkappa_drho = kappa / (2 * rho)
    
    kappaai = kappa *  dii
    kappaai_1 = 1 + kappaai
    
    aux_Xi1 = 1.5 - 2*kappaai_1 + kappaai_1**2/2. 
    aux_Xi1 += np.log(kappaai_1)
    with np.errstate(divide='ignore', invalid='ignore'):
        aux_Xi2 = 3/kappaai**3
        
        daux_Xi1 = (-1 + kappaai + 1/kappaai_1) * dii
        daux_Xi2 = -9 * dii / kappaai**4
        
        Xi = aux_Xi1 * aux_Xi2 
        dXi = aux_Xi1 * daux_Xi2 
        dXi += daux_Xi1 * aux_Xi2 
        
        sum_X = np.dot(Xi, xz2)
        sum_dX = np.dot(dXi, xz2)
        
        aux_a = - betae2_e / (12 * np.pi)
        
        a = kappa * sum_X
        a *= aux_a
        
        da_dkappa = sum_X + kappa * sum_dX
        da_dkappa *= aux_a   
        
        da = da_dkappa * dkappa_drho
    
    return np.nan_to_num(np.array([a, da]))


def d2adh_drho(rho, xz2, er, dii, betae2_e0): 
    betae2_e = betae2_e0 / er
    kappa = np.sqrt(betae2_e * rho * np.sum(xz2))
    dkappa_drho = kappa / (2 * rho)
    d2kappa_drho = - dkappa_drho / (2 * rho)
    
    kappaai = kappa *  dii
    kappaai_1 = 1 + kappaai
    
    aux_Xi1 = 1.5 - 2*kappaai_1 + kappaai_1**2/2. 
    aux_Xi1 += np.log(kappaai_1)
    with np.errstate(divide='ignore', invalid='ignore'):
        aux_Xi2 = 3/kappaai**3
        
        daux_Xi1 = (-1 + kappaai + 1/kappaai_1) * dii
        daux_Xi2 = -9 * dii / kappaai**4
        
        dii2 = dii**2
        d2aux_Xi1 = (1 - 1/kappaai_1**2) * dii2
        d2aux_Xi2 = 36 * dii2 /kappaai**5
        
        Xi = aux_Xi1 * aux_Xi2 
        dXi = aux_Xi1 * daux_Xi2 
        dXi += daux_Xi1 * aux_Xi2 
        
        d2Xi = daux_Xi1 * daux_Xi2 
        d2Xi += aux_Xi1 * d2aux_Xi2 
        d2Xi += d2aux_Xi1 * aux_Xi2 
        d2Xi += daux_Xi1 * daux_Xi2 
        
        aux_a = - betae2_e / (12 * np.pi)
        sum_X = np.dot(Xi, xz2)
        sum_dX = np.dot(dXi, xz2)
        sum_d2X = np.dot(d2Xi, xz2)
        
        a = kappa * sum_X
        a *= aux_a
        
        da_dkappa = sum_X + kappa * sum_dX
        da_dkappa *= aux_a
        
        d2a_dkappa = 2 * sum_dX + kappa * sum_d2X
        d2a_dkappa *= aux_a
        
        
        da = da_dkappa * dkappa_drho
        d2a = d2a_dkappa * dkappa_drho**2
        d2a += da_dkappa * d2kappa_drho
    
    return np.nan_to_num(np.array([a, da, d2a]))

def dadh_dx(rho, derx, z, xz2, er, dii, betae2_e0): 
    betae2_e = betae2_e0 / er
    z2 = z**2
    sum_xz2 = np.sum(xz2)
    kappa = np.sqrt(betae2_e * rho * sum_xz2)
    with np.errstate(divide='ignore', invalid='ignore'):
        dkappa_dx = er * z2 - derx * sum_xz2
        dkappa_dx *= betae2_e0 * rho / (2 * kappa * er**2)
        
        kappaai = kappa *  dii
        kappaai_1 = 1 + kappaai
        
        aux_Xi1 = 1.5 - 2*kappaai_1 + kappaai_1**2/2. 
        aux_Xi1 += np.log(kappaai_1)
    
        aux_Xi2 = 3/kappaai**3
        
        daux_Xi1 = (-1 + kappaai + 1/kappaai_1) * dii
        daux_Xi2 = -9 * dii / kappaai**4
        
        Xi = aux_Xi1 * aux_Xi2 
        dXi = aux_Xi1 * daux_Xi2 
        dXi += daux_Xi1 * aux_Xi2 
        dXi_dx = np.outer(dkappa_dx, dXi)
    
        sum_X = np.dot(Xi, xz2)
        sum_dXdx = dXi_dx@xz2
        sum_dXdx += z2 * Xi
    
        aux_a = - betae2_e / (12 * np.pi)
        a = kappa * sum_X 
        a *= aux_a 
        
        dax = er * sum_dXdx - sum_X * derx
        dax *= kappa
        dax += er * sum_X * dkappa_dx
        dax *= aux_a / er
        
    return np.nan_to_num(a), np.nan_to_num(dax)


def dadh_dxrho(rho, derx, z, xz2, er, dii, betae2_e0): 
    betae2_e = betae2_e0 / er
    z2 = z**2
    sum_xz2 = np.sum(xz2)
    kappa = np.sqrt(betae2_e * rho * sum_xz2)
    dkappa_drho = kappa / (2 * rho)
    dkappa_dx = er * z2 - derx * sum_xz2
    dkappa_dx *= betae2_e0 * rho / (2 * kappa * er**2)
    
    kappaai = kappa *  dii
    kappaai_1 = 1 + kappaai
    
    aux_Xi1 = 1.5 - 2*kappaai_1 + kappaai_1**2/2. 
    aux_Xi1 += np.log(kappaai_1)
    with np.errstate(divide='ignore', invalid='ignore'):
        aux_Xi2 = 3/kappaai**3
        
        daux_Xi1 = (-1 + kappaai + 1/kappaai_1) * dii
        daux_Xi2 = -9 * dii / kappaai**4
        
        Xi = aux_Xi1 * aux_Xi2 
        dXi = aux_Xi1 * daux_Xi2 
        dXi += daux_Xi1 * aux_Xi2 
        dXi_dx = np.outer(dkappa_dx, dXi)   
    
    
        sum_X = np.dot(Xi, xz2)
        sum_dX = np.dot(dXi, xz2)
        sum_dXdx = dXi_dx@xz2
        sum_dXdx += z2 * Xi
    
        aux_a = - betae2_e / (12 * np.pi)
        a = kappa * sum_X 
        a *= aux_a 
        da_dkappa = sum_X + kappa * sum_dX
        da_dkappa *= aux_a   
        da = da_dkappa * dkappa_drho
        
        dax = er * sum_dXdx - sum_X * derx
        dax *= kappa
        dax += er * sum_X * dkappa_dx
        dax *= aux_a / er
        
    return np.nan_to_num(np.array([a, da])), np.nan_to_num(dax)

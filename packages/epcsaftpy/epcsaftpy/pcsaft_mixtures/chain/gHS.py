from __future__ import division, print_function, absolute_import
import numpy as np

def ghs(xhi, aux_dii, aux_dii2):
    xhi0, xhi1, xhi2, xhi3 = xhi
    aux = aux_dii
    aux_2 = aux_dii2
    
    xhi13 = 1 - xhi3
    xhi13_2 = xhi13**2
    xhi13_3 = xhi13_2*xhi13
    
    ghs = 1 / xhi13
    ghs += 3 * aux * xhi2 / xhi13_2
    ghs += 2 * aux_2 * xhi2**2 / xhi13_3
    return ghs


def dghs_dxhi00(xhi, dxhi_dxhi00, aux_dii, aux_dii2):
    xhi0, xhi1, xhi2, xhi3 = xhi
    dxhi0, dxhi1, dxhi2, dxhi3 = dxhi_dxhi00

    aux = aux_dii
    aux_2 = aux_dii2

    xhi2_2 = xhi2**2

    xhi13 = 1 - xhi3
    xhi13_2 = xhi13**2
    xhi13_3 = xhi13_2*xhi13
    xhi13_4 = xhi13_3*xhi13

    ghs = 1 / xhi13
    ghs += 3 * aux * xhi2 / xhi13_2
    ghs += 2 * aux_2 * xhi2_2 / xhi13_3

    dghs = 4. * aux_2 * xhi2 * dxhi2 / xhi13_3
    dghs += 3. * aux * dxhi2 / xhi13_2
    dghs += 6. * aux_2 * xhi2_2 * dxhi3 / xhi13_4
    dghs += 6. * aux * xhi2 * dxhi3 / xhi13_3
    dghs += dxhi3 / xhi13_2
    return np.array([ghs, dghs])


def d2ghs_dxhi00(xhi, dxhi_dxhi00, aux_dii, aux_dii2):
    xhi0, xhi1, xhi2, xhi3 = xhi
    dxhi0, dxhi1, dxhi2, dxhi3 = dxhi_dxhi00
    aux = aux_dii
    aux_2 = aux_dii2

    xhi2_2 = xhi2**2
    dxhi3_2 = dxhi3**2
    xhi13 = 1 - xhi3
    xhi13_2 = xhi13**2
    xhi13_3 = xhi13_2*xhi13
    xhi13_4 = xhi13_3*xhi13
    xhi13_5 = xhi13_4*xhi13

    ghs = 1 / xhi13
    ghs += 3 * aux * xhi2 / xhi13_2
    ghs += 2 * aux_2 * xhi2_2 / xhi13_3

    dghs = 4. * aux_2 * xhi2 * dxhi2 / xhi13_3
    dghs += 3. * aux * dxhi2 / xhi13_2
    dghs += 6. * aux_2 * xhi2_2 * dxhi3 / xhi13_4
    dghs += 6. * aux * xhi2 * dxhi3 / xhi13_3
    dghs += dxhi3 / xhi13_2

    d2ghs = 2 * dxhi3_2 / xhi13_3
    d2ghs += 2 * aux_2 * 2 * dxhi2**2 / xhi13_3
    d2ghs += 2 * aux_2 * 12 * xhi2 * dxhi2 * dxhi3 / xhi13_4
    d2ghs += 2 * aux_2 * 12 * xhi2_2 * dxhi3_2 / xhi13_5
    d2ghs += 3 * aux * 4 * dxhi2 * dxhi3 / xhi13_3
    d2ghs += 3 * aux * 6 * xhi2 * dxhi3_2 / xhi13_4
    return np.array([ghs, dghs, d2ghs])


def dghs_dx(xhi, dxhi_dx, aux_dii, aux_dii2):

    xhi0, xhi1, xhi2, xhi3 = xhi
    dxhi0x, dxhi1x, dxhi2x, dxhi3x = dxhi_dx

    aux = aux_dii
    aux_2 = aux_dii2
    xhi2_2 = xhi2**2

    xhi13 = 1 - xhi3
    xhi13_2 = xhi13**2
    xhi13_3 = xhi13_2*xhi13
    xhi13_4 = xhi13_3*xhi13

    ghs = 1 / xhi13
    ghs += 3 * aux * xhi2 / xhi13_2
    ghs += 2 * aux_2 * xhi2_2 / xhi13_3

    dghsx = 4 * np.multiply.outer(aux_2, dxhi2x) * xhi2 / xhi13_3
    dghsx += 3 * np.multiply.outer(aux, dxhi2x) / xhi13_2
    dghsx += 6 * np.multiply.outer(aux_2, dxhi3x) * xhi2_2 / xhi13_4
    dghsx += 6 * np.multiply.outer(aux, dxhi3x) * xhi2 / xhi13_3
    dghsx += dxhi3x / xhi13_2
    dghsx = np.einsum('ijk->kij',dghsx) 
    #dghsx = dghsx.T

    return ghs, dghsx


def dghs_dxxhi(xhi, dxhi_dxhi00, dxhi_dx, aux_dii, aux_dii2):
    xhi0, xhi1, xhi2, xhi3 = xhi
    dxhi0, dxhi1, dxhi2, dxhi3 = dxhi_dxhi00
    dxhi0x, dxhi1x, dxhi2x, dxhi3x = dxhi_dx

    aux = aux_dii
    aux_2 = aux_dii2
    xhi2_2 = xhi2**2

    xhi13 = 1 - xhi3
    xhi13_2 = xhi13**2
    xhi13_3 = xhi13_2*xhi13
    xhi13_4 = xhi13_3*xhi13
    
    

    ghs = 1 / xhi13
    ghs += 3 * aux * xhi2 / xhi13_2
    ghs += 2 * aux_2 * xhi2_2 / xhi13_3

    dghs = 4. * aux_2 * xhi2 * dxhi2 / xhi13_3
    dghs += 3. * aux * dxhi2 / xhi13_2
    dghs += 6. * aux_2 * xhi2_2 * dxhi3 / xhi13_4
    dghs += 6. * aux * xhi2 * dxhi3 / xhi13_3
    dghs += dxhi3 / xhi13_2   
    
    

    dghsx = 4 * np.multiply.outer(aux_2, dxhi2x) * xhi2 / xhi13_3
    dghsx += 3 * np.multiply.outer(aux, dxhi2x) / xhi13_2
    dghsx += 6 * np.multiply.outer(aux_2, dxhi3x) * xhi2_2 / xhi13_4
    dghsx += 6 * np.multiply.outer(aux, dxhi3x) * xhi2 / xhi13_3
    dghsx += dxhi3x / xhi13_2
   
    dghsx = np.einsum('ijk->kij',dghsx) 

    return np.array([ghs, dghs]), dghsx
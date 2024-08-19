from __future__ import division, print_function, absolute_import
import numpy as np
from .coefficients import aji, bji



def adisp(xhi00, xm, xmi, eta, es3ij, e2s3ij):  
    xm_1 = xm - 1
    xm_2 = xm_1 - 1
    mm1 = xm_1/xm
    mm2 = mm1*xm_2/xm
    ai = aji[0] + aji[1]*mm1 + aji[2]*mm2
    bi = bji[0] + bji[1]*mm1 + bji[2]*mm2
    
    eta2, eta3, eta4, eta5, eta6 = eta**2, eta**3, eta**4, eta**5, eta**6
    etavec = np.array([1, eta, eta2, eta3, eta4, eta5, eta6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
    
    M1 = np.sum(np.outer(xmi,xmi)*es3ij)
    M2 = np.sum(np.outer(xmi,xmi)*e2s3ij)
    
    eta_m1 = 1 - eta
    eta_m2 = 2 - eta
    
    c1n = (20*eta - 27*eta2 + 12*eta3 - 2*eta4)
    c1n *= - xm_1
    c1n /= (eta_m1**2*eta_m2**2)
    c1n += xm * (8*eta - 2*eta2)/eta_m1**4
    c1n += 1
    
    a = - 12 * I1 * M1
    a -= 6 * xm * I2 * M2 / c1n
    a *= xhi00
    return a


def dadisp_dxhi00(xhi00, xm, xmi, eta, deta_dxhi00, es3ij, e2s3ij):
    xm_1 = xm - 1
    xm_2 = xm_1 - 1
    mm1 = xm_1/xm
    mm2 = mm1*xm_2/xm
    ai = aji[0] + aji[1]*mm1 + aji[2]*mm2
    bi = bji[0] + bji[1]*mm1 + bji[2]*mm2
    
    eta2, eta3, eta4, eta5, eta6 = eta**2, eta**3, eta**4, eta**5, eta**6
    etavec = np.array([1, eta, eta2, eta3, eta4, eta5, eta6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
        
    M1 = np.sum(np.outer(xmi,xmi)*es3ij)
    M2 = np.sum(np.outer(xmi,xmi)*e2s3ij)
    
    eta_m1 = eta - 1
    eta_m2 = eta - 2
    
    c1n = (20*eta - 27*eta2 + 12*eta3 - 2*eta4)
    c1n *= - xm_1
    c1n /= (eta_m1**2*eta_m2**2)
    c1n += xm * (8*eta - 2*eta2)/eta_m1**4
    c1n += 1
    
    a = - 12 * I1 * M1
    a -= 6 * xm * I2 * M2 / c1n
    a *= xhi00
    
    ####
    
    detavec = np.array([0, 1, 2*eta, 3*eta2, 4*eta3, 5*eta4, 6*eta5])
    dI1 = np.dot(ai, detavec)
    dI1 *= deta_dxhi00
    dI2 = np.dot(bi, detavec)
    dI2 *= deta_dxhi00
    
    dc1n = 8/eta_m2**3 - 6/eta_m1**3
    dc1n *= - xm_1
    dc1n += xm * (-8 + 4*(-5 + eta)*eta)/eta_m1**5
    dc1n *= deta_dxhi00
    
    da = -xhi00 * I2 * dc1n + c1n * (I2 + xhi00 * dI2)
    da *= - 6 * xm * M2 / c1n**2
    da -= 12 * M1 * (I1 + xhi00 * dI1)
    
    return np.array([a, da])


def d2adisp_dxhi00(xhi00, xm, xmi, eta, deta_dxhi00, es3ij, e2s3ij):
    xm_1 = xm - 1
    xm_2 = xm_1 - 1
    mm1 = xm_1/xm
    mm2 = mm1*xm_2/xm
    ai = aji[0] + aji[1]*mm1 + aji[2]*mm2
    bi = bji[0] + bji[1]*mm1 + bji[2]*mm2
    
    eta2, eta3, eta4, eta5, eta6 = eta**2, eta**3, eta**4, eta**5, eta**6
    etavec = np.array([1, eta, eta2, eta3, eta4, eta5, eta6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
    
    M1 = np.sum(np.outer(xmi,xmi)*es3ij)
    M2 = np.sum(np.outer(xmi,xmi)*e2s3ij)
    
    eta_m1 = eta - 1
    eta_m2 = eta - 2
    
    c1n = (20*eta - 27*eta2 + 12*eta3 - 2*eta4)
    c1n *= - xm_1
    c1n /= (eta_m1**2*eta_m2**2)
    c1n += xm * (8*eta - 2*eta2)/eta_m1**4
    c1n += 1
    
    a = - 12 * I1 * M1
    a -= 6 * xm * I2 * M2 / c1n
    a *= xhi00
    
    ####
    
    detavec = np.array([0, 1, 2*eta, 3*eta2, 4*eta3, 5*eta4, 6*eta5])
    dI1 = np.dot(ai, detavec)
    dI1 *= deta_dxhi00
    dI2 = np.dot(bi, detavec)
    dI2 *= deta_dxhi00
    
    dc1n = 8/eta_m2**3 - 6/eta_m1**3
    dc1n *= - xm_1
    dc1n += xm * (-8 + 4*(-5 + eta)*eta)/eta_m1**5
    dc1n *= deta_dxhi00
    
    da = -xhi00 * I2 * dc1n + c1n * (I2 + xhi00 * dI2)
    da *= - 6 * xm * M2 / c1n**2
    da -= 12 * M1 * (I1 + xhi00 * dI1)
    
    ###
    deta_dxhi00_2 = deta_dxhi00**2
    d2etavec = np.array([0, 0, 2, 6*eta, 12*eta2, 20*eta3, 30*eta4])

    d2I1 = np.dot(ai, d2etavec)
    d2I1 *= deta_dxhi00_2
    d2I2 = np.dot(bi, d2etavec)
    d2I2 *= deta_dxhi00_2
    
    d2c1n = (-2*(-5 + (-6 + eta)*eta)*xm)/eta_m1**6
    d2c1n -= (3/eta_m1**4 - 4/eta_m2**4)*xm_1
    d2c1n *= 6*deta_dxhi00_2
        
    d2a =  -d2c1n*I2*xhi00 - 2*dc1n*(I2 + dI2*xhi00)
    d2a *= c1n
    d2a += c1n**2 * (2*dI2 + d2I2*xhi00)
    d2a += 2*dc1n**2*I2*xhi00
    d2a *= (-6*M2*xm)/c1n**3
    d2a -= 12*M1*(2*dI1 + d2I1*xhi00)
    
    
    
    return np.array([a, da, d2a])
    
def dadisp_dx(x, xhi00, xm, xmi, ms, eta, deta_dx, 
              es3ij, e2s3ij, mes3ij, me2s3ij):
    xm_1 = xm - 1
    xm_2 = xm_1 - 1
    mm1 = xm_1/xm
    mm2 = mm1*xm_2/xm
    ai = aji[0] + aji[1]*mm1 + aji[2]*mm2
    bi = bji[0] + bji[1]*mm1 + bji[2]*mm2
    
    eta2, eta3, eta4, eta5, eta6 = eta**2, eta**3, eta**4, eta**5, eta**6
    etavec = np.array([1, eta, eta2, eta3, eta4, eta5, eta6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
        
    M1 = np.sum(np.outer(xmi,xmi)*es3ij)
    M2 = np.sum(np.outer(xmi,xmi)*e2s3ij)
    
    eta_m1 = eta - 1
    eta_m2 = eta - 2
    
    c1n = (20*eta - 27*eta2 + 12*eta3 - 2*eta4)
    c1n *= - xm_1
    c1n /= (eta_m1**2*eta_m2**2)
    c1n += xm * (8*eta - 2*eta2)/eta_m1**4
    c1n += 1
    
    a = - 12 * I1 * M1
    a -= 6 * xm * I2 * M2 / c1n
    a *= xhi00
    
    ####
     
    aux2_xm = ms/xm**3
    aux1_xm = (-4 + 3*xm)
    aix = np.outer(aji[1]*xm + aji[2]*aux1_xm, aux2_xm)
    bix = np.outer(bji[1]*xm + bji[2]*aux1_xm, aux2_xm)
    
    detavec = np.array([0, 1, 2*eta, 3*eta2, 4*eta3, 5*eta4, 6*eta5])
    dI1x = np.dot(ai, detavec)
    dI1x *= deta_dx
    dI1x += etavec@aix
    dI2x = np.dot(bi, detavec)
    dI2x *= deta_dx
    dI2x += etavec@bix
    
    dc1nx = 8/eta_m2**3 - 6/eta_m1**3
    dc1nx *= - xm_1
    dc1nx += xm * (-8 + 4*(-5 + eta)*eta)/eta_m1**5
    dc1nx *= deta_dx
    dc1nx += - ms*(20*eta - 27*eta2 + 12*eta3 - 2*eta4)/(eta_m1**2*eta_m2**2)
    dc1nx += ms * (8*eta - 2*eta2)/eta_m1**4
    
    dM1x = 2*mes3ij@x
    dM2x = 2*me2s3ij@x
    
    dax = ms * I2 * M2 + dM2x * I2 * xm 
    dax += dI2x * M2 * xm - (dc1nx * I2 * M2 * xm) / c1n
    dax *= -6/c1n
    dax -= 12 * (I1 * dM1x + dI1x * M1)
    dax *= xhi00
    
    return a, dax


def dadisp_dxxhi(x, xhi00, xm, xmi, ms, eta, deta_dx, 
                 deta_dxhi00, es3ij, e2s3ij, mes3ij, me2s3ij):
    xm_1 = xm - 1
    xm_2 = xm_1 - 1
    mm1 = xm_1/xm
    mm2 = mm1*xm_2/xm
    ai = aji[0] + aji[1]*mm1 + aji[2]*mm2
    bi = bji[0] + bji[1]*mm1 + bji[2]*mm2
    
    eta2, eta3, eta4, eta5, eta6 = eta**2, eta**3, eta**4, eta**5, eta**6
    etavec = np.array([1, eta, eta2, eta3, eta4, eta5, eta6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
        
    M1 = np.sum(np.outer(xmi,xmi)*es3ij)
    M2 = np.sum(np.outer(xmi,xmi)*e2s3ij)
    
    eta_m1 = eta - 1
    eta_m2 = eta - 2
    
    c1n = (20*eta - 27*eta2 + 12*eta3 - 2*eta4)
    c1n *= - xm_1
    c1n /= (eta_m1**2*eta_m2**2)
    c1n += xm * (8*eta - 2*eta2)/eta_m1**4
    c1n += 1
    
    a = - 12 * I1 * M1
    a -= 6 * xm * I2 * M2 / c1n
    a *= xhi00
    
    ####
    
    detavec = np.array([0, 1, 2*eta, 3*eta2, 4*eta3, 5*eta4, 6*eta5])
    dI1 = np.dot(ai, detavec)
    dI1 *= deta_dxhi00
    dI2 = np.dot(bi, detavec)
    dI2 *= deta_dxhi00
    
    dc1n = 8/eta_m2**3 - 6/eta_m1**3
    dc1n *= - xm_1
    dc1n += xm * (-8 + 4*(-5 + eta)*eta)/eta_m1**5
    dc1n *= deta_dxhi00
    
    da = -xhi00 * I2 * dc1n + c1n * (I2 + xhi00 * dI2)
    da *= - 6 * xm * M2 / c1n**2
    da -= 12 * M1 * (I1 + xhi00 * dI1)
    
    ####
     
    aux2_xm = ms/xm**3
    aux1_xm = (-4 + 3*xm)
    aix = np.outer(aji[1]*xm + aji[2]*aux1_xm, aux2_xm)
    bix = np.outer(bji[1]*xm + bji[2]*aux1_xm, aux2_xm)
    
    detavec = np.array([0, 1, 2*eta, 3*eta2, 4*eta3, 5*eta4, 6*eta5])
    dI1x = np.dot(ai, detavec)
    dI1x *= deta_dx
    dI1x += etavec@aix
    dI2x = np.dot(bi, detavec)
    dI2x *= deta_dx
    dI2x += etavec@bix
    
    dc1nx = 8/eta_m2**3 - 6/eta_m1**3
    dc1nx *= - xm_1
    dc1nx += xm * (-8 + 4*(-5 + eta)*eta)/eta_m1**5
    dc1nx *= deta_dx
    dc1nx += - ms*(20*eta - 27*eta2 + 12*eta3 - 2*eta4)/(eta_m1**2*eta_m2**2)
    dc1nx += ms * (8*eta - 2*eta2)/eta_m1**4
    
    dM1x = 2*mes3ij@x
    dM2x = 2*me2s3ij@x
    
    dax = ms * I2 * M2 + dM2x * I2 * xm 
    dax += dI2x * M2 * xm - (dc1nx * I2 * M2 * xm) / c1n
    dax *= -6/c1n
    dax -= 12 * (I1 * dM1x + dI1x * M1)
    dax *= xhi00
    
    return np.array([a, da]), dax
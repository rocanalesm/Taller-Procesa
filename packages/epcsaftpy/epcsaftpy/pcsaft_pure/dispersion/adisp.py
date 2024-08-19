from __future__ import division, print_function, absolute_import
import numpy as np


def adisp(ms, eta, dia3, ai, bi, m2e2s3, mes3):
    c1 = 1/(1 + ((20*eta - 27*eta**2 + 12*eta**3 - 2*eta**4)*(1 - ms))/((1 - eta)**2*(2 - eta)**2) + ((8*eta - 2*eta**2)*ms)/(1 - eta)**4)
    
    etavec = np.array([1, eta, eta**2, eta**3, eta**4, eta**5, eta**6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
    
    a = -6*eta*(c1*I2*m2e2s3 + 2*I1*mes3)/dia3
    return a


def dadisp_deta(ms, eta, dia3, ai, bi, m2e2s3, mes3):
    c1 = 1/(1 + ((20*eta - 27*eta**2 + 12*eta**3 - 2*eta**4)*(1 - ms))/((1 - eta)**2*(2 - eta)**2) + ((8*eta - 2*eta**2)*ms)/(1 - eta)**4)
    dc1 = ((-2*(-1 + eta)**2*(20 + eta*(-24 + eta*(6 + eta))) - 2*(12 + (-1 + eta)*eta*(-96 + eta*(90 + (-25 + eta)*eta)))*ms)/((-2 + eta)**3*(-1 + eta)**5))*c1**2
    
    etavec = np.array([1, eta, eta**2, eta**3, eta**4, eta**5, eta**6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
    
    detavec = np.array([0, 1, 2*eta, 3*eta**2, 4*eta**3, 5*eta**4, 6*eta**5])
    dI1 = np.dot(ai, detavec)
    dI2 = np.dot(bi, detavec)
    
    factor = -6*eta/dia3
    a = factor*(c1*I2*m2e2s3 + 2*I1*mes3)
    da = factor*(c1*dI2*m2e2s3 + dc1*I2*m2e2s3 + 2*dI1*mes3) + a/eta
    return np.array([a, da])


def d2adisp_deta(ms, eta, dia3, ai, bi, m2e2s3, mes3):
    c1 = 1/(1 + ((20*eta - 27*eta**2 + 12*eta**3 - 2*eta**4)*(1 - ms))/((1 - eta)**2*(2 - eta)**2) + ((8*eta - 2*eta**2)*ms)/(1 - eta)**4)
    dc1 = ((-2*(-1 + eta)**2*(20 + eta*(-24 + eta*(6 + eta))) - 2*(12 + (-1 + eta)*eta*(-96 + eta*(90 + (-25 + eta)*eta)))*ms)/((-2 + eta)**3*(-1 + eta)**5))*c1**2
    d2c1 = (2*(-((-1 + eta)**6*(-1072 + eta*(3936 + eta*(-6456 + eta*(6200 + eta*(-3744 + eta*(1368 + eta*(-250 + 3*eta*(2 + eta))))))))) + ms*(-1 + eta)**4*(528 + eta*(8208 + eta*(-39252 + eta*(77200 + eta*(-88460 + eta*(65596 + eta*(-32137 + eta*(10004 + eta*(-1756 + 3*eta*(40 + eta)))))))))) + ms**2*(-1 + eta)**2*(576 + eta*(5040 + eta*(8172 + eta*(-95280 + eta*(233112 + eta*(-311268 + eta*(269557 + eta*(-160048 + eta*(65442 + eta*(-17812 + eta*(2971 + 6*(-42 + eta)*eta)))))))))))))/(4 + eta*(-(eta*(26 + eta*(-42 + eta*(27 + (-8 + eta)*eta)))) + ms*(12 + eta*(27 + eta*(-70 + eta*(51 + 2*(-8 + eta)*eta))))))**3

    etavec = np.array([1, eta, eta**2, eta**3, eta**4, eta**5, eta**6])
    I1 = np.dot(ai, etavec)
    I2 = np.dot(bi, etavec)
    
    detavec = np.array([0, 1, 2*eta, 3*eta**2, 4*eta**3, 5*eta**4, 6*eta**5])
    dI1 = np.dot(ai, detavec)
    dI2 = np.dot(bi, detavec)
    
    d2etavec = np.array([0, 0, 2, 6*eta, 12*eta**2, 20*eta**3, 30*eta**4])
    d2I1 = np.dot(ai, d2etavec)
    d2I2 = np.dot(bi, d2etavec)
    
    factor = -6*eta/dia3
    a = factor*(c1*I2*m2e2s3 + 2*I1*mes3)
    da = factor*(c1*dI2*m2e2s3 + dc1*I2*m2e2s3 + 2*dI1*mes3) + a/eta
    d2a = -(2*(6*c1*dI2*m2e2s3 + 6*dc1*I2*m2e2s3 + 12*dI1*mes3) + (6*(c1*d2I2 + 2*dc1*dI2 + d2c1*I2)*m2e2s3 + 12*d2I1*mes3)*eta)/dia3

    return np.array([a, da, d2a])

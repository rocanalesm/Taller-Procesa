from __future__ import division, print_function, absolute_import
import numpy as np



def ahs(eta):
    a = ((4 - 3*eta)*eta)/(eta - 1)**2
    return a


def dahs_deta(eta):
    a = ((4 - 3*eta)*eta)/(eta - 1)**2
    da = (2*(-2 + eta))/(eta - 1)**3
    return np.array([a, da])


def d2ahs_deta(eta):
    a = ((4 - 3*eta)*eta)/(eta - 1)**2
    da = (2*(-2 + eta))/(eta - 1)**3
    d2a = (10 - 4*eta)/(eta - 1)**4
    return np.array([a, da, d2a])
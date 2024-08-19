from __future__ import division, print_function, absolute_import
import numpy as np


def ghs(eta):
    g = (-2 + eta)/(2.*(eta - 1)**3)
    return g


def dghs_deta(eta):
    g = (-2 + eta)/(2.*(eta - 1)**3)
    dg = (5 - 2*eta)/(2.*(eta - 1)**4)
    return np.array([g, dg])


def d2ghs_deta(eta):
    g = (-2 + eta)/(2.*(eta - 1)**3)
    dg = (5 - 2*eta)/(2.*(eta - 1)**4)
    d2g = (3*(-3 + eta))/(eta - 1)**5
    return np.array([g, dg, d2g])
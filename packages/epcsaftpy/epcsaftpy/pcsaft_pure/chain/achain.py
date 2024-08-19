from __future__ import division, print_function, absolute_import
import numpy as np


def achain(ms, ghs):
    a = (1 - ms)*np.log(ghs)
    return a


def dachain_deta(ms, dghs):
    dghs_ghs = dghs[1]/dghs[0]
    a = (1 - ms)*np.log(dghs[0])
    da = (1 - ms)*dghs_ghs
    return np.array([a, da])


def d2achain_deta(ms, d2ghs):
    dghs_ghs = d2ghs[1]/d2ghs[0]
    d2ghs_ghs = d2ghs[2]/d2ghs[0]
    a = (1 - ms)*np.log(d2ghs[0])
    da = (1 - ms)*dghs_ghs
    d2a = (1 - ms)*(d2ghs_ghs - dghs_ghs**2)
    return np.array([a, da, d2a])

from __future__ import division, print_function, absolute_import
from .pcsaft_mixtures.pcsaftmix import pcsaft_mix
from .pcsaft_pure.pcsaft import pcsaft_pure



def pcsaft(mix_or_component, compute_critical=True):
    '''
    Returns PC-SAFT EoS object.

    Parameters
    ----------
    mix_or_component : object
        :class:`epcsaftpy.mixture` or :class:`epcsaftpy.component` object
    compute_critical: bool
        If True the critical point of the fluid will attempt to be computed
        (it might fail for some fluids).
    Returns
    -------
    eos : object
        PC-SAFT EoS object
    '''
    nc = mix_or_component.nc
    if nc == 1:
        eos = pcsaft_pure(mix_or_component, compute_critical)
    else:
        eos = pcsaft_mix(mix_or_component, compute_critical)
    return eos




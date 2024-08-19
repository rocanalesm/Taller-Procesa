
# epcsaftpy

epcsaftpy provide you the PC-SAFT EoS to perform phase equilibria and interfacial properties using [sgtpy python package](https://github.com/gustavochm/sgtpy). This module includes PC-SAFT EoS with the following contributions for the Helmholtz energy:

```math
   a = a^{id} + a^{hc} + a^{disp} + a^{assoc} + a^{polar} + a^{born} + a^{DH} 
```

where

- $a^{id}$: ideal contribution
- $a^{hc}$: hard chain contribution
- $a^{disp}$: dispersion contribution
- $a^{assoc}$: association contribution
- $a^{polar}$: polar contribution (using [Jog & Chapman term](https://doi.org/10.1080/00268979909482832))
- $a^{born}$: Born contribution
- $a^{DH}$: Debye-Hückel contribution

**New in v0.0.09**: json file as a input for component definition is extended for electrolyte components.

**New in v0.0.08**: A bug in the electrolyte implementation has been resolved, deleting the dispersion terms exclusively in cation-cation and anion-anion interactions when the ions have a spherical shape ($ms = 1$). Also, the adjustment of temperature-dependent diameters ($d = \sigma * (1 - 0.12)$) applies exclusively to ions with a spherical shape. This improves the modeling of ionic liquids with electrolyte terms.

**New in v0.0.07**: json file as a input for component definition.

**New in v0.0.06**: Free volume theory and Helmholtz scaling to obtain the viscosity.

**New in v0.0.05**: NET-GP implementation to consider the thermodynamics of glassy mixtures.

**New in v0.0.04**: Numba implementation to improve the performance in the association contribution.

**New in v0.0.03**: Entropy scaling to obtain the viscosity.


## Installation Prerequisites

Necessary

- numpy
- scipy
- numba


Optional (only for examples)

- sgtpy
- matplotlib
- xlsxwriter
- pandas
- python-ternary

## Installation


To get the git beta version, run:
```console
$ git clone https://github.com/estebancea/epcsaftpy
$ cd epcsaftpy
$ pip install .
```
    
## Getting Started

First, components are defined with their molecular parameters, then a mixture can be created with them.

```python
from epcsaftpy import component, pcsaft
Water = component('Water', ms = 1.2046817736, sigma = [2.7927, 10.11, -0.01775, -1.417, -0.01146], eps = 353.9449,
                 kappaAB = 0.045090, eAB = 2425.6714, sites = [0, 1, 1], Mw = 18.01528)
FA = component('Furfuryl Alcohol', ms = 4.361081, sigma = 3.004829 , eps = 218.33885, 
               kappaAB = 0.14622, eAB = 1834.334, sites = [0, 1, 2], Mw = 98.1014)
mix = Water + FA
# adding the binary parameters
mix.set_kijsaft(i = 0, j = 1, kij0 = -0.01)
eos = pcsaft(mix)
```
The eos object can be used to compute phase equilibria using [sgtpy module](https://github.com/gustavochm/sgtpy).

```python
import numpy as np
from sgtpy.equilibrium import bubblePy
# computing bubble point
T = 298.15 # K
x = np.array([0.3, 0.7])
# initial guesses for vapor composition and pressure
y0 = 1.*x
P0 = 8000. # Pa
sol = bubblePy(y0, P0, x, T, eos, full_output = True)
```
Finally, the equilibria results can be used to model the interfacial behavior of the mixture [sgtpy module](https://github.com/gustavochm/sgtpy).

For more examples, please have a look at the Jupyter Notebook files
located in the *examples* folder of the sources

## Documentation

For more details about how to use the epcsaftpy module, see [this documentation](https://github.com/estebancea/epcsaftpy/tree/main/docs).

## License information
See `LICENSE.md` for information on the terms & conditions for usage
of this software, and a DISCLAIMER OF ALL WARRANTIES.

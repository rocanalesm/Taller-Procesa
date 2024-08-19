
## Component object

Creates an object with pure component parameters and then saved as class attributes.

```python
comp1 = epcsaftpy.component(name='None', ms=1., sigma=0., eps=0., eAB=0.,
                            kappaAB=0, sites=[0, 0, 0], mupol=0, xpol=0., 
                            z = 0, er = 0, Mw=1., cii=0.):
```
### Parameters
`name` : *str*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Name of the component

`ms` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Chain lenght

`sigma` : *float* or *list*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Segment diameter, input in [Ángstrom], stored in [m].\
&emsp; &emsp;&emsp;&emsp;&emsp;For thermal sigma, it should be a list like: $\left[\sigma_0,\ t_1,\ t_2,\ t_3,\ t_4\right]$ where $\sigma = \sigma_0 + t_1 \cdot \exp \left(t_2 \cdot T \right) + t_3 \cdot \exp \left( t_4 \cdot T \right)$.

`eps` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Dispersion energy, input in [K], stored in [J].

`eAB` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Association Energy, input in [K], stores in [J].

`kappaAB` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Association volume.

`sites` : *list*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Association sites [Bivalent, Positive, Negative].

`mupol` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;dipolar moment [Debye].

`xpol` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;fraction of dipolar segment on a chain molecule (Jog and Chapman polar term).

`z` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Ions valence.

`er` : *float* or *list*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Dielectric constant.\
&emsp; &emsp;&emsp;&emsp;&emsp;For thermal dielectric constant, it should be a list like: $\left[\varepsilon_{r0},\ \varepsilon_{r1}\right]$ where $\varepsilon_{r} = \varepsilon_{r0} + \varepsilon_{r1} \cdot T$.

`Mw` : *float*,\
&emsp; &emsp;&emsp;&emsp;&emsp;Molar weight [g $\cdot$ mol<sup>-1</sup>].

`cii` : *list*,\
&emsp; &emsp;&emsp;&emsp;&emsp;polynomial coefficient for influence parameter in SGT for sgtpy module [J $\cdot$ m<sup>5</sup> $\cdot$ mol<sup>-2</sup>].

## Mixture object
The parameters for each component and the mixing rules are created as an object to be saved as class attributes. A mixture object is created adding component objects as follow:

```python
mix = comp1 + comp2 + comp3 + comp4
```

The parameters for the mixing rules are added using the following methods for a mixture object.
### Methods
```python
mix.set_kijsaft(i, j, kij0, kij1 = 0., kij2 = 0., kij3 = 0.)
```
&emsp;&emsp;Method that sets the kij correction for dispersion energy between component i and component j.

```math
    \epsilon_{ij} = (1 - k_{ij}) \sqrt{\epsilon_i \epsilon_j}
``` 
&emsp;&emsp;kij correction is computed as follows:

```math
    k_{ij} = k_{ij,0} + k_{ij,1} \cdot T +  k_{ij,2} \cdot T^2 + k_{ij,3} / T
```


```python
mix.set_kepsijsaft(i, j, kepsij0, kepsij1 = 0., kepsij2 = 0., kepsij3 = 0.)
```
&emsp;&emsp;Method that sets the kepsij correction for association energy between component i and component j.

```math
   \epsilon_{AiBj} = (1 - keps_{ij}) \dfrac{\epsilon_{Ai} + \epsilon_{Bj}}{2}
``` 
&emsp;&emsp;kepsij correction is computed as follows:

```math
    keps_{ij} = keps_{ij,0} + keps_{ij,1} \cdot T +  keps_{ij,2} \cdot T^2 + keps_{ij,3} / T
```

```python
mix.set_lijsaft(i, j, lij0, lij1 = 0., lij2 = 0., lij3 = 0.)
```
&emsp;&emsp;Method that sets the lij correction for segment diameter between component i and component j.

```math
   \sigma_{ij} = (1 - l_{ij}) \dfrac{\sigma_{i} + \sigma_{j}}{2}
``` 
&emsp;&emsp;lij correction is computed as follows:

```math
    l_{ij} = l_{ij,0} + l_{ij,1} \cdot T +  l_{ij,2} \cdot T^2 + l_{ij,3} / T
```

[Return to previous page](https://github.com/estebancea/epcsaftpy/tree/main/docs)

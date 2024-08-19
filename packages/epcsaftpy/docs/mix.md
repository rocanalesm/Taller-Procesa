
## EoS for mixtures

Pure component PC-SAFT EoS Object
```python
eos = epcsaftpy.pcsaft(mixture, compute_critical = True)
```
This object have implemeted methods for phase equilibrium as for interfacial properties calculations based in [sgtpy module](https://github.com/gustavochm/sgtpy).

### Parameters

`mixture` : *object*, \
&emsp;&emsp;&emsp;Mixture created with mixture class.

`compute_critical` : *bool*, \
&emsp;&emsp;&emsp;If True the critical point for each pure component will attempt to be computed (it might fail for some fluids).

### Attributes

`eos.nc`: *integrer*, \
&emsp;&emsp;&emsp;Number of component in the mixture.

`eos.Mw`: *array_like*, \
&emsp;&emsp;&emsp;Molar weight for each component [g $\cdot$ mol<sup>-1</sup>].

`eos.sigma`: *array_like*, \
&emsp;&emsp;&emsp;Segment diameter for each component [m].

`eos.eps`: *array_like*, \
&emsp;&emsp;&emsp;Dispersion energy for each component [J].

`eos.eAB`: *array_like*, \
&emsp;&emsp;&emsp;Association energy for each component [J].

`eos.kappaAB`: *array_like*, \
&emsp;&emsp;&emsp;Association volume for each component.

`eos.sigmaij`: *array_like*, \
&emsp;&emsp;&emsp;Segment diameter matrix [m].

```math
\sigma_{ij} = \dfrac{\sigma_i + \sigma_j}{2} \left(1 - l_{ij}\right)
```

`eos.epsij`: *array_like*, \
&emsp;&emsp;&emsp;Dispersion energy matrix [J]. 

```math
\varepsilon_{ij} = \sqrt{\varepsilon_{i} \varepsilon_{j}} \left(1 - k_{ij}\right)
```
`eos.eABij`: *array_like*, \
&emsp;&emsp;&emsp;Association energy matrix [J].

```math 
\varepsilon_{AiBj} = \dfrac{\varepsilon_{Ai} + \varepsilon_{Bj}}{2} \left(1 - keps_{ij}\right)
```
`eos.kappaABij`: *array_like*, \
&emsp;&emsp;&emsp;Association volume.

```math 
\kappa_{AiBj} = \sqrt{\kappa_{Ai} \kappa_{Bj}} \left( \dfrac{\sqrt{\sigma_i  \sigma_j}}{\sigma_{ij}}\right)^3
```

`eos.sitesmix`: *list*, \
&emsp;&emsp;&emsp;Triplet of number of association sites [Bivalents, Positives, Negatives] for each component.

`eos.mupol`: *array_like*, \
&emsp;&emsp;&emsp;Dipolar moment [Debye].

`eos.xpol`: *array_like*, \
&emsp;&emsp;&emsp;Fraction of dipolar segment sites. 

`eos.cii`: *array_like*, \
&emsp;&emsp;&emsp;Influence factor for SGT [J $\cdot$ m<sup>5</sup> $\cdot$ mol<sup>-2</sup>]. 

`eos.cij`: *array_like*, \
&emsp;&emsp;&emsp;Cross influence parameter matrix for SGT [J $\cdot$ m<sup>5</sup> $\cdot$ mol<sup>-2</sup>].

`eos.beta`: *array_like*, \
&emsp;&emsp;&emsp;Correction to cross influence parameter matrix.

**Warning**: *Some attributes are temperature dependence, so they are only defined after a temperature calculation and overwriting after a new one.*      
        

### Methods 

```python
eos.density(x, T, P, state, rho0 = None, Xass0 = None)
```
&emsp;&emsp;Computes the density of the of the mixture at given composition, temperature, pressure and aggregation state (liquid or vapor). 

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`rho0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute density root [mol $\cdot$ m<sup>-3</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp;&emsp;&emsp;**Return**

&emsp;&emsp;&emsp;`rho` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar density [mol $\cdot$ m<sup>-3</sup>].
</details>

```python
eos.pressure(x, rho, T, Xass0 = None)
```
&emsp;&emsp;Computes the pressure at given composition, density and temperature. 

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`rho` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp; Molar density [mol $\cdot$ m<sup>-3</sup>].

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp;&emsp;&emsp;**Return**

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].
</details>

```python
eos.logfugef(x, T, P, state, v0 = None, Xass0 = None)
```
&emsp;&emsp;Computes the logarithmus of the effective fugacity coefficient of the components in the mixture at given composition, temperature, pressure and aggregation state (liquid or vapor). 

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>

&emsp; &emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`v0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute volume root [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp;&emsp;&emsp;**Returns**

&emsp;&emsp;&emsp;`lnphi` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Logarithmus of the effective fugacity coefficient for each component.

&emsp;&emsp;&emsp;`v` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar volume [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].
</details>

```python
eos.get_lnphi_pure(T, P, state)
```
&emsp;&emsp;Computes the logarithm of the pure component's fugacity coefficient at given temperature, pressure and aggregation state (liquid or vapor).

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp; &emsp;&emsp;**Parameters**

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp; &emsp;&emsp;`P` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp; &emsp;&emsp;`state` : *string*, \
&emsp; &emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp; &emsp;&emsp;**Returns**

&emsp; &emsp;&emsp;`lnphi_pure` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Logarithmus of the effective fugacity coefficient for each component.

</details>

```python
eos.get_lngamma(x, T, P, v0 = None, Xass0 = None, lnphi_pure = None)
```
&emsp;&emsp;Computes computes the activity coefficient of the mixture at given molar composition, temperature, pressure and aggregation state (liquid or vapor).

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp; &emsp;&emsp;**Parameters**

&emsp; &emsp;&emsp;`x` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp; &emsp;&emsp;`P` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp; &emsp;&emsp;`v0` : *float, optional*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess to compute volumen root [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp; &emsp;&emsp;`lnphi_pure` : *array_like, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp; Logarithm of the pure components's fugacity coefficient. Computed if not provided.

&emsp; &emsp;&emsp;**Returns**

&emsp; &emsp;&emsp;`lnphi_pure` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Logarithmus of the effective fugacity coefficient for each component.
</details>

```python
eos.muad(rhoi, T, Xass0 = None)
```
&emsp;&emsp;Computes the dimentionless chemical potential at given density vector and temperature.

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
              
&emsp; &emsp;&emsp;**Parameters**
        
&emsp; &emsp;&emsp;`rhoi` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Molar density vector $\rho_i = x_i \rho$ [mol $\cdot$ m<sup>-3</sup>].
        
&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 
        
&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp; &emsp;&emsp;**Return**
        
&emsp; &emsp;&emsp;`muad` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Dimentionless chemical potencial. 
</details>

```python
eos.dmuad(rhoi, T, Xass0 = None)
```
&emsp;&emsp;Computes the chemical potential and its numerical derivative at given density vector and temperature.

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
              
&emsp; &emsp;&emsp;**Parameters**
        
&emsp; &emsp;&emsp;`rhoi` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Molar density vector $\rho_i = x_i \rho$ [mol $\cdot$ m<sup>-3</sup>].
        
&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 
        
&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp; &emsp;&emsp;**Return**
        
&emsp; &emsp;&emsp;`muad` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Dimentionless chemical potencial. 
        
        
&emsp; &emsp;&emsp;`dmuad` : *array_like*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Derivavites of the dimentionless chemical potencial respect to rhoi [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].
</details>


```python
eos.EntropyR(x, T, P, state, v0 = None, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Computes the residual entropy of the mixture at given composition, temperature, pressure and aggregation state (liquid or vapor).

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`v0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute volume root [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp;&emsp;&emsp;`T_step` : *float, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy [K]. 
        
&emsp;&emsp;&emsp;**Return**

&emsp;&emsp;&emsp;`Sr` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Residual entropy [J/mol K].
        
</details>

```python
eos.EnthalpyR(x, T, P, state, v0=None, Xass0=None, T_step = 0.1)
```
&emsp;&emsp;Computes the residual enthalpy of the mixture at given composition, temperature, pressure and aggregation state (liquid or vapor).

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`v0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute volume root [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp;&emsp;&emsp;`T_step` : *float, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy [K]. 
        
&emsp;&emsp;&emsp;**Return**
        
&emsp;&emsp;&emsp;`Hr` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Residual enthalpy [J/mol].
</details>

```python
eos.CvR(x, rho, T, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Computes the residual isochoric heat capacity of the mixture at given composition, density and temperature.

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.
        
&emsp;&emsp;&emsp;`rho` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar density [mol $\cdot$ m<sup>-3</sup>].

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp;&emsp;&emsp;`T_step` : *float, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy [K]. 
        
&emsp;&emsp;&emsp;**Return**
        
&emsp;&emsp;&emsp;`Cv` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Residual isochoric heat capacity [J/mol K].
        
</details>


```python
eos.CpR(x, T, P, state, v0 = None, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Computes the residual heat capacity of the mixture at given composition, temperature, pressure and aggregation state (liquid or vapor). 

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`v0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute volume root [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp;&emsp;&emsp;`T_step` : *float, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy [K]. 
        
&emsp;&emsp;&emsp;**Return**
        
&emsp;&emsp;&emsp;`Cp` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Residual heat capacity [J/mol K].
        
</details>

```python
eos.speed_sound(x, T, P, state, v0 = None, Xass0 = None, T_step = 0.1,
                CvId = 3*R/2, CpId = 5*R/2)
```
&emsp;&emsp;Computes the speed of sound of the mixture at given composition, temperature, pressure and aggregation state (liquid or vapor). *This calculation requires that the molar weight of the fluids has been set in the component function.* 

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
        
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`x` : *array_like*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar fraction array.

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K]. 

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`v0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute volume root [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
        
&emsp; &emsp;&emsp;`T_step` : *float, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy [K]. 
        
&emsp; &emsp;&emsp;`CvId` : *float, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Ideal gas isochoric heat capacity, set to 3R/2 by default [J/mol K].
        
&emsp; &emsp;&emsp;`CpId` : *float, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Ideal gas heat capacity, set to 3R/2 by default [J/mol K].
        
&emsp;&emsp;&emsp;**Return**
        
&emsp; &emsp;&emsp;`w` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;Speed of sound [m/s].

</details>

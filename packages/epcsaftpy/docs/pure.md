## EoS for pure component

Pure component PC-SAFT EoS Object
```python
eos = epcsaftpy.pcsaft(pure, compute_critical = True)
```
This object have implemeted methods for phase equilibrium as for interfacial properties calculations based in [sgtpy module](https://github.com/gustavochm/sgtpy).


### Parameters  

`pure` : *object*, \
&emsp;&emsp;&emsp;Pure component created with component class.

`compute_critical` : *bool*, \
&emsp;&emsp;&emsp;If True the critical point of the fluid will attempt to be computed (it might fail for some fluids).


### Attributes 

`eos.ms` : *float*, \
&emsp;&emsp;&emsp;Number of chain segments. 

`eos.sigma` : *float*, \
&emsp;&emsp;&emsp;Segment diameter parameter [m]. 

`eos.eps` : *float*, \
&emsp;&emsp;&emsp;Dispersion energy [J]. 

`eos.eABij` : *float*, \
&emsp;&emsp;&emsp;Association energy [J]. 

`eos.kappaABij` : *float*, \
&emsp;&emsp;&emsp;Association volume. 

`eos.sites`: *list*, \
&emsp;&emsp;&emsp;Triplet of number of association sites [Bivalents, Positives, Negatives].

`eos.mupol`: *float*, \
&emsp;&emsp;&emsp;Dipolar moment [Debye]. 

`eos.xpol`: *float*, \
&emsp;&emsp;&emsp;Fraction of dipolar segment on a chain molecule.

`eos.cii`: *list or array*, \
&emsp;&emsp;&emsp;Influence factor for SGT [J $\cdot$ m<sup>5</sup> $\cdot$ mol<sup>-2</sup>]. 


### Methods <br />
```python
eos.density(T, P, state, rho0 = None, Xass0 = None)
```
&emsp;&emsp;A method that computes fluid density at a given temperature and pressure.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp;&emsp;&emsp;**Parameters**
 
&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K].

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Pressure [Pa].

&emsp;&emsp;&emsp;`state` : *string*, \
&emsp;&emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase.

&emsp;&emsp;&emsp;`rho0` : *float, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess to compute density root [mol $\cdot$ m<sup>-3</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp;&emsp;&emsp;**Returns**

&emsp;&emsp;&emsp;`rho` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Molar density [mol $\cdot$ m<sup>-3</sup>].
</details>

```python
eos.psat(T, P0 = None, v0 = [None, None], Xass0 = [None, None], full_output = False)
```
&emsp;&emsp;A method that computes saturation pressure at a fixed T.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`T` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Absolute temperature [K].

&emsp;&emsp;&emsp;`P0` : *float, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;initial value to find saturation pressure [Pa].

&emsp;&emsp;&emsp;`v0` : *list, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;initial guess for liquid and vapor phase, respectively [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp;&emsp;&emsp;`full_output` : *bool, optional*, \
&emsp;&emsp;&emsp;&emsp;&emsp;If you want to output all the results of the calculations.

&emsp;&emsp;&emsp;**Returns**

&emsp;&emsp;&emsp;`psat` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;saturation pressure [Pa].

&emsp;&emsp;&emsp;`vl`: *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;liquid saturation volume [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`vv`: *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;vapor saturation volume [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].
</details>

```python
eos.tsat(P, T0 = None, Tbounds = None, v0 = [None, None], Xass0 = [None, None], full_output = False)
```
&emsp;&emsp;A method that computes saturation temperature at a given pressure.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp;&emsp;&emsp;**Parameters**

&emsp;&emsp;&emsp;`P` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;absolute pressure [Pa].

&emsp;&emsp;&emsp;`T0` : *float, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;Temperature to start iterations [K].

&emsp;&emsp;&emsp;`Tbounds` : *tuple, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;(Tmin, Tmax) Temperature interval to start iterations [K].

&emsp;&emsp;&emsp;`v0` : *list, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;initial guess for liquid and vapor phase, respectively [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`Xass0` : *array, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp;&emsp;&emsp;`full_output` : *bool, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;Whether to outputs or not all the calculation info.

&emsp;&emsp;&emsp;**Returns**

&emsp;&emsp;&emsp;`tsat` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Saturation temperature [K]. 

&emsp;&emsp;&emsp;`vl` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Liquid saturation volume [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].

&emsp;&emsp;&emsp;`vv` : *float*, \
&emsp;&emsp;&emsp;&emsp;&emsp;Vapor saturation volume [m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].
</details>

```python
eos.get_critical(Tc0 = None, rhoc0 = None, method = 'hybr', full_output = False, overwrite = False)
```
&emsp;&emsp;A method that solves the critical coordinate of the fluid using good initial guesses for the critical temperature and density as inputs:

<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp;&emsp;&emsp;**Parameters**
 
&emsp;&emsp;&emsp;`Tc0` : *float, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;initial guess for critical temperature [K].

&emsp;&emsp;&emsp;`rhoc` : *float, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;initial guess for critical density [mol $\cdot$ m<sup>-3</sup>].

&emsp;&emsp;&emsp;`method` : *string, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;SciPy; root method to solve critical coordinate.

&emsp;&emsp;&emsp;`full_output` : *bool, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;whether to outputs or not all the calculation info.

&emsp;&emsp;&emsp;`overwrite` : *bool, optional* \
&emsp;&emsp;&emsp;&emsp;&emsp;wheter to overwrite already computed critical points.

&emsp;&emsp;&emsp;**Returns**

&emsp;&emsp;&emsp;`Tc` : *float*, \
&emsp;&emsp;&emsp;Critical temperature [K].

&emsp;&emsp;&emsp;`Pc` : *float*, \
&emsp;&emsp;&emsp;Critical pressure [Pa]

&emsp;&emsp;&emsp;`rhoc` : *float*, \
&emsp;&emsp;&emsp;Critical density [mol $\cdot$ m<sup>-3</sup>]
</details>


```python
eos.afcn(rho, T, Xass0 = None)
```
&emsp; &emsp;Method that computes the total Helmholtz free energy of the fluid.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters**

 
&emsp; &emsp;&emsp;`rho` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;molecular density [molecules $\cdot$ m<sup>-3</sup>].

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K].

&emsp; &emsp;&emsp;`Xass0`: *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp; &emsp;&emsp;**Returns**

&emsp; &emsp;&emsp;`a`: *float*, \
&emsp; &emsp;&emsp;Helmholtz free energy [J$\cdot$ mol<sup>-1</sup>].
</details>

```python
eos.dafcn_drho(rho, T, Xass0 = None)
```
&emsp;&emsp; Method that computes the total Helmholtz free energy of the fluid and its first density derivative. 
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters**

&emsp; &emsp;&emsp;`rho` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;molecular density [molecules $\cdot$ m<sup>-3</sup>].

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K].

&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.

&emsp; &emsp;&emsp;**Returns**

&emsp; &emsp;&emsp;a: float, <br />
&emsp; &emsp;&emsp;Helmholtz free energy and its derivative  [J $\cdot$ mol<sup>-1</sup>, J m<sup>3</sup> $\cdot$ mol<sup>-1</sup>].
</details>

```python
eos.pressure(rho, T, Xass0 = None)
```
&emsp;&emsp;Method that computes the pressure at given density [mol $\cdot$ m<sup>-3</sup>] and temperature [K]
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters**
 
&emsp; &emsp;&emsp;`rho` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;molecular density [molecules $\cdot$ m<sup>-3</sup>].

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K].

&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
 
&emsp; &emsp;&emsp;**Returns**
 
&emsp; &emsp;&emsp;`P` : *float*, \
&emsp; &emsp;&emsp;Pressure [Pa].
 
</details>

```python
eos.dP_drho(rho, T, Xass0 = None)
```
&emsp;&emsp;Method that computes the pressure and its density derivative at given density [mol $\cdot$ m<sup>-3</sup>] and temperature [K].
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;rho: float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;density [mol $\cdot$ m<sup>-3</sup>] <br />
&emsp; &emsp;&emsp;T : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;P : float <br />
&emsp; &emsp;&emsp;Pressure [Pa] <br />
&emsp; &emsp;&emsp;dP: float <br />
&emsp; &emsp;&emsp;derivate of pressure respect density [Pa m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />
</details>

```python
eos.logfug(T, P, state, v0 = None, Xass0 = None)
```
&emsp;&emsp;Method that computes the fugacity coefficient at given composition, temperature and pressure.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;T : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K]<br />
&emsp; &emsp;&emsp; P : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;pressure [Pa] <br />
&emsp; &emsp;&emsp;state : string <br />
&emsp; &emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase <br />
&emsp; &emsp;&emsp;v0: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;initial guess for volume root [m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;logfug : float <br />
&emsp; &emsp;&emsp;fugacity coefficient <br />
&emsp; &emsp;&emsp;v : float <br />
&emsp; &emsp;&emsp;computed volume of the phase [m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />

</details>

```python
eos.a0ad(rho, T, Xass0 = None)
```
&emsp;&emsp;Method that computes the adimenstional Helmholtz density energy at given density and temperature.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters**
 
&emsp; &emsp;&emsp;`rho` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;molecular density [molecules $\cdot$ m<sup>-3</sup>].

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K].

&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
 
&emsp; &emsp;&emsp;**Returns**
 
&emsp; &emsp;&emsp;`a0ad`:  *float*, \
&emsp; &emsp;&emsp;Helmholtz density energy [J $\cdot$ m<sup>-3</sup>].
 
</details>


```python
eos.muad(rho, T, Xass0 = None)
```
&emsp;&emsp;Method that computes the adimenstional chemical potential at given density and temperature.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters**
 
&emsp; &emsp;&emsp;`rho` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;molecular density [molecules $\cdot$ m<sup>-3</sup>].

&emsp; &emsp;&emsp;`T` : *float*, \
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K].

&emsp; &emsp;&emsp;`Xass0` : *array, optional* \
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites.
 
&emsp; &emsp;&emsp;**Returns**
 
&emsp; &emsp;&emsp;`muad` : *float*, \
&emsp; &emsp;&emsp;Chemical potential [J $\cdot$ mol<sup>-1</sup>]:
</details>


```python
eos.dOm(rho, T, mu, Psat, Xass0 = None)
```
&emsp;&emsp;Method that computes the adimenstional Thermodynamic Grand potential at given density and temperature.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;rho : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;density [mol $\cdot$ m<sup>-3</sup>] <br />
&emsp; &emsp;&emsp;T : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;mu : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;adimentional chemical potential at equilibrium <br />
&emsp; &emsp;&emsp;Psat : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;adimentional pressure [Pa] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;Out: float <br />
&emsp; &emsp;&emsp;Thermodynamic Grand potential [Pa] <br />
</details>


```python
eos.ci(T)
```
&emsp;&emsp;Method that evaluates the polynomial for the influence parameters used in the SGT theory for surface tension calculations.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;T : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;ci: float <br />
&emsp; &emsp;&emsp;influence parameters [J m<sup>5</sup> mol<sup>-2</sup>] <br />
</details>


```python
eos.EntropyR(T, P, state, v0 = None, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Method that computes the residual entropy at given temperature and pressure.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;T : float<br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;P : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;pressure [Pa] <br />
&emsp; &emsp;&emsp;state : string <br />
&emsp; &emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase <br />
&emsp; &emsp;&emsp;v0: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;initial guess for volume root [m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br />
&emsp; &emsp;&emsp;T_step: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;Sr : float <br />
&emsp; &emsp;&emsp;residual entropy [J $\cdot$ mol<sup>-1</sup> K<sup>-1</sup>] <br />
</details>



```python
eos.EnthalpyR(T, P, state, v0 = None, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Method that computes the residual entropy at given temperature and pressure.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;T : float<br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;P : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;pressure [Pa] <br />
&emsp; &emsp;&emsp;state : string <br />
&emsp; &emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase <br />
&emsp; &emsp;&emsp;v0: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;initial guess for volume root [m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br />
&emsp; &emsp;&emsp;T_step: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;Hr : float <br />
&emsp; &emsp;&emsp;residual enthalpy [J $\cdot$ mol<sup>-1</sup>] <br />
</details>


```python
eos.CvR(rho, T, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Method that computes the residual isochoric heat capacity at given density and temperature.
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;rho : float<br />
&emsp; &emsp;&emsp;&emsp;&emsp;density [mol $\cdot$ m<sup>-3</sup>] <br />
&emsp; &emsp;&emsp;T : float<br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br />
&emsp; &emsp;&emsp;T_step: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;Cv: float <br />
&emsp; &emsp;&emsp;isochoric heat capacity [J $\cdot$ mol<sup>-1</sup> K<sup>-1</sup> ] <br />
</details>

```python
eos.CpR(T, P, state, v0 = None, Xass0 = None, T_step = 0.1)
```
&emsp;&emsp;Method that computes the residual isochoric heat capacity at given density and temperature.
 <details closed>
 &emsp;&emsp;&emsp;<summary>  
Show details:
 </summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;T : float<br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;P : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;pressure [Pa] <br />
&emsp; &emsp;&emsp;state : string <br />
&emsp; &emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase <br />
&emsp; &emsp;&emsp;v0: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;initial guess for volume root [m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br />
&emsp; &emsp;&emsp;T_step: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;Cp: float <br />
&emsp; &emsp;&emsp;isochoric heat capacity [J $\cdot$ mol<sup>-1</sup> K<sup>-1</sup> ] <br />
</details>

```python
eos.speed_sound(T, P, state, v0 = None, Xass0 = None, T_step = 0.1,
                CvId = 3 * R / 2, CpId = 5 * R / 2)
```
&emsp;&emsp;Method that computes the speed of sound at given temperature and pressure. This calculation requires that the molar weight of the fluid has been set in the component function. By default the ideal gas Cv and Cp are set to 3R/2 and 5R/2, the user can supply better values if available.
  
<details closed>
&emsp;&emsp;&emsp;<summary>  
Show details:
</summary>
 
&emsp; &emsp;&emsp;**Parameters** <br />
&emsp; &emsp;&emsp;T : float<br />
&emsp; &emsp;&emsp;&emsp;&emsp;absolute temperature [K] <br />
&emsp; &emsp;&emsp;P : float <br />
&emsp; &emsp;&emsp;&emsp;&emsp;pressure [Pa] <br />
&emsp; &emsp;&emsp;state : string <br />
&emsp; &emsp;&emsp;&emsp;&emsp;'L' for liquid phase and 'V' for vapor phase <br />
&emsp; &emsp;&emsp;v0: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;initial guess for volume root [m<sup>3</sup> $\cdot$ m<sup>-1</sup>] <br />
&emsp; &emsp;&emsp;Xass0: array, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Initial guess for the calculation of fraction of non-bonded sites <br />
&emsp; &emsp;&emsp;T_step: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Step to compute the numerical temperature derivates of Helmholtz free energy <br /> 
&emsp; &emsp;&emsp;CvId: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Ideal gas isochoric heat capacity, set to 3R/2 by default [J $\cdot$ mol<sup>-1</sup> K<sup>-1</sup> ] <br /> 
&emsp; &emsp;&emsp;CpId: float, optional <br />
&emsp; &emsp;&emsp;&emsp;&emsp;Ideal gas heat capacity, set to 3R/2 by default [J $\cdot$ mol<sup>-1</sup> K<sup>-1</sup> ] <br /> 
&emsp; &emsp;&emsp;**Returns** <br />
&emsp; &emsp;&emsp;w: float <br />
&emsp; &emsp;&emsp;speed of sound [m $\cdot$ s<sup>-1</sup>] <br />
</details>


[Return to previous page](https://github.com/estebancea/epcsaftpy/tree/main/docs)


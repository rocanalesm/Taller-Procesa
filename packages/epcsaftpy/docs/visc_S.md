## Theory: entropy scaling

The Chapman−Enskog viscosity of substance $i$ ( $\eta_{CE,i}$ ) is defined by:

```math
\eta_{CE,i} = \dfrac{5}{16} \dfrac{\sqrt{Mw_i k_B T / (N_A \pi )}}{\sigma_i^2 \Omega_i^{(2,2)}}
```
where the collision integral $\Omega_i^{(2,2)}$ is calculated using the reduced temperature ( $T_{ad,i}=k_B T / \varepsilon_i$ ) by:

```math
    \Omega_i^{(2,2)} = 1.16145 T_{ad,i}^{- 0.14874} + 0.52487 \exp(-0.7732 T_{ad,i}) + 2.16178 \exp(-2.43787 T_{ad,i}) -6.435\cdot 10^{-4}  T_{ad,i}^{0.14874}  \sin(18.0323  T_{ad,i}^{-0.7683} - 7.27371)
```

This viscosity is used as a reference to formulate a dimensionless form of the viscosity $\eta_{i}^{*}=\eta_{i}/\eta_{CE,i}$, which is fitted for each pure component to a polynomial that is a function of the reduced residual entropy as follow:

```math
\ln \eta_i = A_i + B_i s_i^* + C_i s_i^{*2}+ D_i s_i^{*3}
```

where $s_i^*$ is the dimensionless residual entropy:

```math
s_i^* = \dfrac{s_{res}\left(x_i=1, \rho, T \right)}{N_A k_B m_i}
```

**Obs.**: Some parameters in the literature were fitted to a different representation of the Chapman−Enskog viscosity (considering the segment number). In this case, only one parameter has to be changed in this way: $A_i = A_{i,old} + \ln \left(1/\sqrt{m_i} \right)$.

### Mixtures
Dimensionless viscosity for mixtures:
```math
\ln \eta_i = \sum_{i = 1}^n x_i A_i + \sum_{i = 1}^n \dfrac{x_i m_i}{\bar{m}} B_i s^* + \sum_{i = 1}^n \dfrac{x_i m_i}{\bar{m}} C_i s^{*2}+ \sum_{i = 1}^n \dfrac{x_i m_i}{\bar{m}} D_i s^{*3}
```

with

```math
\bar{m} = \sum_{i = 1}^n x_i m_i
```
```math
s^* = \dfrac{s_{res}\left(x, \rho, T \right)}{N_A k_B \bar{m}}
```

The Chapman−Enskog viscosity for mixtures:
```math
\eta_{CE,mix} = \sum_{i = 1}^n \dfrac{x_i \eta_{CE,i}}{\sum_{j = 1}^{n} x_j \phi_{ij}}
```

where

```math
\phi_{ij}=\dfrac{\left(1 + \left(\eta_{CE,i} / \eta_{CE,j}\right)^{1/2} \left(Mw_{j} / Mw_{i}\right)^{1/4}  \right)^2}{\left(8\left(1 + Mw_{i} / Mw_{j}\right)\right)^{1/2}}
```

### References
[Ind. Eng. Chem. Res. 2018, 57, 4095−4114](https://pubs.acs.org/doi/full/10.1021/acs.iecr.7b04871)

[Return to previous page](https://github.com/estebancea/epcsaftpy/tree/main/docs)


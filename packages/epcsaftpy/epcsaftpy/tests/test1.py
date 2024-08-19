import numpy as np
from epcsaftpy import component, pcsaft 

print("\n#################################################################")
print("|        Checking the mixtures in each function                  |")     
print("#################################################################")

Na = 6.022142e23
tol = 1e-10
comp1 = component('comp1', ms=1.96720036, sigma=4.54762477, eps=377.60127994, xpol = 0.9, mupol = 2.5, 
                  eAB = 2500., kappaAB = 0.0024, sites = [0, 1, 1], er = 50. ) 
comp2 = component('comp2', ms=3.96720036, sigma=3.54762477, eps=233.60127994, xpol = 0.9, mupol = 1.5, 
                  eAB = 1530., kappaAB = 0.041, sites = [1, 2, 1], er = 80. ) 
comp3 = component('comp3', ms=2.96720036, sigma=2.54762477, eps=133.60127994, xpol = 0.5, mupol = 3.5, 
                  eAB = 2004., kappaAB = 0.002, sites = [2, 1, 5], z = 1, er = 8.)
comp4 = component('comp4', ms=12.96720036, sigma=2.54762477, eps=133.60127994, xpol = 0.33, mupol = 2.445, 
                  eAB = 4234., kappaAB = 0.032, sites = [2, 10, 5], z = 2, er = 8.)
mix = comp1 + comp2 + comp3 + comp4
mix.set_kijsaft(i=0, j=1, kij0=0.003, kij1=0.0004, kij2=0.00005, kij3=0.006)
mix.set_kijsaft(i=0, j=2, kij0=0.02)
mix.set_kijsaft(i=1, j=2, kij0=-0.03)
mix.set_lijsaft(i=0, j=1, lij0=-0.01)
mix.set_lijsaft(i=0, j=2, lij0=0.003)
mix.set_lijsaft(i=1, j=2, lij0=0.006, lij1=0.0007, lij2=0.00002, lij3=0.002)
mix.set_kepsijsaft(i=0, j=1, kepsij0=-0.01)
mix.set_kepsijsaft(i=0, j=2, kepsij0=0.043)
mix.set_kepsijsaft(i=1, j=3, kepsij0=0.046, kepsij1=0.0007, kepsij2=0.00002, kepsij3=0.002)
eos = pcsaft(mix)


def check_fun(fun1, fun2, name_fun1, name_fun2):
    space = 13
    space1 = (space - len(name_fun1))*" "
    space2 = (space - len(name_fun2))*" "
    if np.min(np.abs(fun2)) != 0:
        check = np.average(abs((fun1 - fun2)/fun2))
    else:
        check = np.average(abs((fun1 - fun2)))
    if check==0:
        log10 = -99
    else:
        log10 =  int(round(np.log10(check), 0))
    #print("Checking   " + name_fun1 + " == " + name_fun2 
    #      + f"       : log10(E) = {log10} |", check < tol)
    print("Checking {}{:>} == {:<}{}  | log10(E) = {} | {} |".format(space1, name_fun1, name_fun2, 
                                                                     space2, log10, check < tol))
    
def check_contribution(Contribution, da, dax):
    d2a_drhoF, da_dxF = Contribution
    aF, da_drhoF, d2a_drhoF = d2a_drhoF
    a, da_drho, d2a_drho = da
    da_dx = dax
    print("_________________________________________________________________")
    check_fun(aF, a, "aF", "a")
    check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
    check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
    check_fun(da_dxF, da_dx, "da_dxF", "da_x")

def checking(T, rho, x):
    rhom = Na * rho
    # Check all the function are the same
    aresF = eos.ares(x, rhom, T)
    dares_drhoF = eos.dares_drho(x, rhom, T)
    d2ares_drhoF = eos.d2ares_drho(x, rhom, T)
    dares_dxF1, dares_dxF2 = eos.dares_dx(x, rhom, T)
    dares_dxrhoF1, dares_dxrhoF2 = eos.dares_dxrho(x, rhom, T)



    check_fun(dares_drhoF[0], aresF, "dares_drho", "ares")
    check_fun(d2ares_drhoF[0], aresF, "d2ares_drho", "ares")
    check_fun(dares_dxF1, aresF, "dares_dx", "ares")
    check_fun(dares_dxrhoF1[0], aresF, "d2ares_dxrho", "ares")
    check_fun(d2ares_drhoF[0:2], dares_drhoF, "d2ares_drho", "dares_drho")
    check_fun(dares_dxF2, dares_dxrhoF2, "dares_dxrho", "dares_drho")
    check_fun(dares_dxrhoF2, dares_dxF2, "dares_dxrho", "dares_dx")

 

print("_________________________________________________________________")
T = 250.
rho = 1/1.1
x = np.array([0.2, 0.1, 0.5, 0.2])
checking(T, rho, x)
print("_________________________________________________________________")
T = 600.
rho = 1/20.
x = np.array([0.5, 0.2, 0.2, 0.1])
checking(T, rho, x)
print("_________________________________________________________________")

print("\n#################################################################")
print("|        Checking the mixtures in the pure component              |")     
print("#################################################################")
## Checking the pure component
# Componets parameters
comp1 = component('comp1', ms=1.96720036, sigma=4.54762477, eps=377.60127994, 
                  xpol = 0.8, mupol = 1.5, 
                  eAB = 2500., kappaAB = 0.04, sites = [0, 1, 1])
comp2 = component('comp2', ms=3.96720036, sigma=[2.7927, 10.11 , -0.01775, -1.417 , -0.01146], eps=233.60127994, 
                  xpol = 0.9, mupol = 2.5)
comp3 = component('comp3', ms=2.96720036, sigma=2.54762477, eps=133.60127994)

mix = comp1 + comp2 + comp3
mix.set_kijsaft(i=0, j=1, kij0=0.003, kij1=0.0004, kij2=0.00005, kij3=0.006)
mix.set_kijsaft(i=0, j=2, kij0=0.02)
mix.set_kijsaft(i=1, j=2, kij0=-0.03)
mix.set_lijsaft(i=0, j=1, lij0=-0.01)
mix.set_lijsaft(i=0, j=2, lij0=0.003)
mix.set_lijsaft(i=1, j=2, lij0=0.006, lij1=3.0e-5, lij2=2.0e-6, lij3=2.0e-5)
eos = pcsaft(mix)
eos1 = pcsaft(comp1)
eos2 = pcsaft(comp2)
eos3 = pcsaft(comp3)

rho = 0.9
rhom = Na * rho
T = 300.0

print("\n Checking the mixture in the pure component \n  * No Associative component")
x = np.array([0., 0., 1.])
aM, da_drhoM, d2a_drhoM = eos.d2ares_drho(x, rhom, T)
[aP, da_drhoP, d2a_drhoP], _ = eos3.d2ares_drho(rhom, T)
print("_________________________________________________________________")
check_fun(aM, aP, "aM", "aP")
check_fun(da_drhoM, da_drhoP, "da_drhoM", "da_drhoP")
check_fun(d2a_drhoM, d2a_drhoP, "d2a_drhoM", "d2a_drhoP")

print("\n Checking the mixture in the pure component \n  * Associative component")
x = np.array([1., 0., 0.])
aM, da_drhoM, d2a_drhoM = eos.d2ares_drho(x, rhom, T)
[aP, da_drhoP, d2a_drhoP], _ = eos1.d2ares_drho(rhom, T)
print("_________________________________________________________________")
check_fun(aM, aP, "aM", "aP")
check_fun(da_drhoM, da_drhoP, "da_drhoM", "da_drhoP")
check_fun(d2a_drhoM, d2a_drhoP, "d2a_drhoM", "d2a_drhoP")


print("\n Checking the mixture in the pure component \n  * Polar component")
x = np.array([0., 1., 0.])
aM, da_drhoM, d2a_drhoM = eos.d2ares_drho(x, rhom, T)
[aP, da_drhoP, d2a_drhoP], _ = eos2.d2ares_drho(rhom, T)
print("_________________________________________________________________")
check_fun(aM, aP, "aM", "aP")
check_fun(da_drhoM, da_drhoP, "da_drhoM", "da_drhoP")
check_fun(d2a_drhoM, d2a_drhoP, "d2a_drhoM", "d2a_drhoP")



print("\n#################################################################")
print("|          Ternary Mixture ( Different Points )                 |")     
print("#################################################################")
# Componets parameters
comp1 = component('comp1', ms=1.96720036, sigma=4.54762477, eps=377.60127994, 
                  xpol = 0.8, mupol = 1.5, 
                  eAB = 2500., kappaAB = 0.04, sites = [0, 1, 1])
comp2 = component('comp2', ms=3.96720036, sigma=3.54762477, eps=233.60127994, 
                  xpol = 0.9, mupol = 2.5)
comp3 = component('comp3', ms=2.96720036, sigma=2.54762477, eps=133.60127994, 
                  xpol = 0.5, mupol = 3.5)



## Checking with Wolfram Mathematica (POINT 1)
mix = comp1 + comp2 + comp3
eos = pcsaft(mix)

rho = 1/(1.2)
rhom = Na * rho
T = 200.0
x = np.array([0.3, 0.2, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.045624098902866055, -8.184530539568543e-26, 2.9745125797211347e-50
da_dx = np.array([-0.25705184586980445, -0.013622515194335636, -0.014280771984983029])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")


## Checking with Wolfram Mathematica (POINT 2)
mix = comp1 + comp2 + comp3
mix.set_kijsaft(i=0, j=1, kij0=0.02)
mix.set_kijsaft(i=0, j=2, kij0=0.01)
mix.set_kijsaft(i=1, j=2, kij0=-0.03)
eos = pcsaft(mix)

rho = 1/(1.2)
rhom = Na * rho
T = 200.0
x = np.array([0.3, 0.2, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.045614985068071275, -8.182714728563173e-26, 2.9745115672275425e-50
da_dx = np.array([-0.2569884635704341, -0.013597660547847924, -0.014284763580901822])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")


## Checking with Wolfram Mathematica (POINT 3)
mix = comp1 + comp2 + comp3
mix.set_kijsaft(i=0, j=1, kij0=0.003, kij1=0.0004, kij2=0.00005, kij3=0.006)
mix.set_kijsaft(i=0, j=2, kij0=0.01)
mix.set_kijsaft(i=1, j=2, kij0=-0.03)
eos = pcsaft(mix)

rho = 1/(1.2)
rhom = Na * rho
T = 200.0
x = np.array([0.3, 0.2, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.04524838234206607, -8.109658880629436e-26, 2.974531036150099e-50
da_dx = np.array([-0.25585163710466735,-0.011936474947798594, -0.014413297609127941])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")

## Checking with Wolfram Mathematica (POINT 4)
mix = comp1 + comp2 + comp3
mix.set_lijsaft(i=0, j=1, lij0=0.02)
mix.set_lijsaft(i=0, j=2, lij0=0.01)
mix.set_lijsaft(i=1, j=2, lij0=-0.03)
eos = pcsaft(mix)

rho = 1/(1.2)
rhom = Na * rho
T = 200.0
x = np.array([0.3, 0.2, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.04560822385909381, -8.181367556885784e-26, 2.974511175722919e-50
da_dx = np.array([-0.25692910094450644, -0.013589188830565261, -0.014294435038156722])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")

## Checking with Wolfram Mathematica (POINT 5)
mix = comp1 + comp2 + comp3
mix.set_kijsaft(i=0, j=1, kij0=0.003, kij1=0.0004, kij2=0.00005, kij3=0.006)
mix.set_kijsaft(i=0, j=2, kij0=0.02)
mix.set_kijsaft(i=1, j=2, kij0=-0.03)
mix.set_lijsaft(i=0, j=1, lij0=-0.01)
mix.set_lijsaft(i=0, j=2, lij0=0.003)
mix.set_lijsaft(i=1, j=2, lij0=0.006, lij1=0.0007, lij2=0.00002, lij3=0.002)
eos = pcsaft(mix)

rho = 1/(1.2)
rhom = Na * rho
T = 200.0
x = np.array([0.3, 0.2, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.04501751607505954, -8.063657744390406e-26, 2.9745216754985556e-50
da_dx = np.array([-0.2558149194483608, -0.010794684602093565, -0.013927071878044305])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")

## Checking with Wolfram Mathematica (POINT 6)
comp1 = component('comp1', ms=1.96720036, sigma=4.54762477, eps=377.60127994, 
                  xpol = 0.8, mupol = 1.5, 
                  eAB = 2500., kappaAB = 0.04, sites = [0, 1, 1])
comp2 = component('comp2', ms=3.96720036, sigma=[2.7927, 10.11 , -0.01775, -1.417 , -0.01146], eps=233.60127994, 
                  xpol = 0.9, mupol = 2.5)
comp3 = component('comp3', ms=2.96720036, sigma=2.54762477, eps=133.60127994, 
                  xpol = 0.5, mupol = 3.5)




mix = comp1 + comp2 + comp3
mix.set_kijsaft(i=0, j=1, kij0=0.003, kij1=0.0004, kij2=0.00005, kij3=0.006)
mix.set_kijsaft(i=0, j=2, kij0=0.02)
mix.set_kijsaft(i=1, j=2, kij0=-0.03)
mix.set_lijsaft(i=0, j=1, lij0=-0.01)
mix.set_lijsaft(i=0, j=2, lij0=0.003)
mix.set_lijsaft(i=1, j=2, lij0=0.006, lij1=3e-5, lij2=2e-6, lij3=2e-5)
eos = pcsaft(mix)


rho = 1/(1.2)
rhom = Na * rho
T = 200.0
x = np.array([0.3, 0.2, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.04577683001083403, -8.214898220120929e-26, 2.974734155633182e-50
da_dx = np.array([-0.25583739110224196, -0.015062408732569251, -0.015160687844158645])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")



## Checking with Wolfram Mathematica (POINT 7)

T = 250.0
x = np.array([0.7, 0.2, 0.1])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.0223237161473961443, -4.34094570341571640e-26, 4.08850253013830908e-51
da_dx = np.array([-0.0609288538228760918, -0.00904133004989607185, -0.00673737646973781613])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")

## Checking with Wolfram Mathematica (POINT 8)

rho = 0.9
rhom = Na * rho
T = 150.0
x = np.array([0.25, 0.25, 0.5])
# Wolfram Mathematica Point
a, da_drho, d2a_drho = -0.460240046951552184, -3.50141802991918388e-25, 5.07679905580106073e-49
da_dx = np.array([-2.51918197208977190, -0.0254557940384418820, -0.0285170191681427959])
# Python Point
d2a_drhoF = eos.d2ares_drho(x, rhom, T)
_, da_dxF = eos.dares_dx(x, rhom, T)
aF, da_drhoF, d2a_drhoF = d2a_drhoF

print("_________________________________________________________________")
check_fun(aF, a, "aF", "a")
check_fun(da_drhoF, da_drho, "da_drhoF", "da_drho")
check_fun(d2a_drhoF, d2a_drho, "d2a_drhoF", "d2a_drho")
check_fun(da_dxF, da_dx, "da_dxF", "da_dx")
print("_________________________________________________________________")


print("\n#################################################################")
print("|          Ternary Mixture ( + A_assoc + A_polar )               |")     
print("#################################################################")

rho = 0.9
rhom = Na * rho
T = 150.0
x = np.array([0.25, 0.25, 0.5])

Mono, Chain, Disp,  AA, Polar, _, _ = eos.contributions(x, rhom, T)

print("\n Checking all the value with Wolfram Mathematica \n  * Monomers contribution")
da = 0.000322881384766006668, 5.95751287255215182e-28, 7.83779550297340936e-56
dax = np.array([0.000761277234524941225, 0.000885891841166293391, 0.000467964024422461253])
check_contribution(Mono, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Chain contribution")
da = -0.000134629207462323187, -2.48400828233796709e-28, -1.53054379763763565e-56
dax = np.array([-0.000291762431823127321, -0.000392081236384062297, -0.000196599491718527460])
check_contribution(Chain, da, dax)



print("\n Checking all the value with Wolfram Mathematica \n  * Dispersion contribution")
da = -0.00157323472018657332, -2.90240684004992414e-27, 1.02914935729128733e-54
dax = np.array([-0.00521780946909759964, -0.00390232360356286122, -0.00253971994144350710])
check_contribution(Disp, da, dax)


print("\n Checking all the value with Wolfram Mathematica \n  * Association contribution")
da = -0.449299023402536524, -3.29957835973934668e-25, 5.07669972414204373e-49
dax = np.array([-2.51249591304939024, -0.0000212707170158228221, -8.98452646361745915e-6])
check_contribution(AA, da, dax)


print("\n Checking all the value with Wolfram Mathematica \n  * Polar contribution")
da = -0.00955604100613275179, -1.76289106369552254e-26, 8.84094402728366016e-54
dax = np.array([-0.00193776437398611489, -0.0220260103226454299, -0.0262396792329396035])
check_contribution(Polar, da, dax)




#####################
print("\n#################################################################")
print("|          Ternary Mixture ( + A_DH + A_Born )                  |")     
print("#################################################################")
 
# Componets parameters
comp1 = component('Methanol', ms = 1.5255208343918, sigma = 3.2300, eps = 188.904644, 
                  eAB = 2899.49055, kappaAB = 0.0351760892, sites = [0, 1, 1], 
                  er = [90.087, -0.1921], Mw = 32.042)
comp2 = component('Li', ms = 1.0, sigma = 2.8449294, eps = 360.00, z = 1., er = 8., Mw = 6.941)
comp3 = component('Br', ms = 1.0, sigma = 3.0707333, eps = 190.00, z = -1., er = 8., Mw = 79.904000)

mix = comp1 + comp2 + comp3
eos = pcsaft(mix)
rho = 0.9
rhom = Na * rho
T = 298.15
x = np.array([0.6, 0.2, 0.2])
Mono, Chain, Disp,  Asso, Polar, DH, Born = eos.contributions(x, rhom, T)

print("\n Checking all the value with Wolfram Mathematica \n  * Monomers contribution")
da = 0.00005359362846436879, 9.888381827955725e-29, 4.678100576734754e-57
dax = np.array([0.0001345648839078139, 0.00006262277818304839, 0.00006962229024711415])
check_contribution(Mono, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Chain contribution")
da = -8.324895711002956e-6, -1.535987625676060e-29, -3.153786249487104e-58
dax = np.array([-0.00002476727155120895, -4.072686879695217e-6, -4.874687185981040e-6])
check_contribution(Chain, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Dispersion contribution")
da = -0.0001364604466082001, -2.517739331512263e-28, 5.396418445505396e-57
dax = np.array([-0.0003015494268920141, -0.0002006763805632452, -0.0001553868748224704])
check_contribution(Disp, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Dispersion contribution")
da = -0.0001364604466082001, -2.517739331512263e-28, 5.396418445505396e-57
dax = np.array([-0.0003015494268920141, -0.0002006763805632452, -0.0001553868748224704])
check_contribution(Disp, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Debye-Huckel contribution")
da = -0.02625071196284421, -2.383757238604602e-26, 2.267875608024107e-50
dax = np.array([0.05615648435189368, -0.08427396079446364, -0.08419549226113002])
check_contribution(DH, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Born contribution")
da = -82.48740614787209, 0., 0.
dax = np.array([-5.402983084768301, -215.4072859387275, -199.6643515003987])
check_contribution(Born, da, dax)
print("#################################################################")



#####################
print("\n#################################################################")
print("|          Binary Mixture ( + A_DH + A_Born )                  |")     
print("#################################################################")
 
# Componets parameters
cation = component('C2min', ms = 1.4872, sigma = 3.5926, eps = 206.49, z = 1., er = 1., Mw = 111.1675)
anion = component('BF4', ms = 3.8227, sigma = 3.5088, eps = 496.12, z = -1., er = 1., Mw = 86.805)

mix = cation + anion
eos = pcsaft(mix)

rho = 5. # mol/m3
rhom = Na * rho
T = 298.15
x = np.array([0.5, 0.5])
Mono, Chain, Disp,  Asso, Polar, DH, Born = eos.contributions(x, rhom, T)

print("\n Checking all the value with Wolfram Mathematica \n  * Monomers contribution")
da = 0.0019311350078041972, 6.414906724730165e-28, 9.685022488981736e-56
dax = np.array([0.0021855079073576098, 0.00553991002955878])
check_contribution(Mono, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Chain contribution")
da = -0.0007517028319619825, -2.4967131243633662e-28, -1.6594182408283351e-56
dax = np.array([-0.0006502495318250988, -0.0023567122289168516])
check_contribution(Chain, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Dispersion contribution")
da = -0.011244079480124756, -3.732563067628561e-27, 1.1174953715912088e-54
dax = np.array([-0.008942330008232777, -0.03991176973100235])
check_contribution(Disp, da, dax)

print("\n Checking all the value with Wolfram Mathematica \n  * Debye-Huckel contribution")
da = -19.75718437654341, -2.396738650490602e-24, 6.087255193924686e-49
dax = np.array([0.024695217112025603, -0.024695217112014944])
check_contribution(DH, da, dax)
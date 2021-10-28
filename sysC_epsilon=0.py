from math import pi

import numpy as np
from astropy.constants import G, M_sun, c
from scipy.integrate import quad


def prod_11(f):
    return 51370.1092733744/(f**(7/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_12(f):
    return 0

def prod_13(f):
    return 0

def prod_14(f):
    return 0

def prod_15(f):
    return 0

def prod_16(f):
    return 0

def prod_17(f):
    return 0

def prod_22(f):
    return 41.9347830803056*pi**2/(f**(1/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_23(f):
    return -1467.7174078107*pi/(f**(4/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_24(f):
    return -74652884.9824338*(0.024105978576213*pi**(2/3)*f**(2/3) - 0.00315231420649041*pi**2*f + 1)/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_25(f):
    return -19907435.9953157*(0.0147197441049309*pi**(2/3)*f**(2/3) - 0.0106390604469051*pi**2*f)/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_26(f):
    return 88248.5062071077*pi**(1/3)/(f**2*(f**2/2450 + 2 + 24010000/f**4))

def prod_27(f):
    return -17423.3052960052*pi**(2/3)/(f**(5/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_33(f):
    return 51370.1092733744/(f**(7/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_34(f):
    return 2612850974.38518*(0.024105978576213*pi**(2/3)*f**(2/3) - 0.00315231420649041*pi**2*f + 1)/(pi**(5/3)*f**4*(f**2/2450 + 2 + 24010000/f**4))

def prod_35(f):
    return 696760259.836049*(0.0147197441049309*pi**(2/3)*f**(2/3) - 0.0106390604469051*pi**2*f)/(pi**(5/3)*f**4*(f**2/2450 + 2 + 24010000/f**4))

def prod_36(f):
    return -3088697.71724877/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_37(f):
    return 609815.685360182/(pi**(1/3)*f**(8/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_44(f):
    return 132898105745009.0*(0.024105978576213*pi**(2/3)*f**(2/3) - 0.00315231420649041*pi**2*f + 1)**2/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_45(f):
    return 35439494865335.8*(0.0147197441049309*pi**(2/3)*f**(2/3) - 0.0106390604469051*pi**2*f)*(0.024105978576213*pi**(2/3)*f**(2/3) - 0.00315231420649041*pi**2*f + 1)/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_46(f):
    return -157101220033.372*(0.024105978576213*pi**(2/3)*f**(2/3) - 0.00315231420649041*pi**2*f + 1)/(pi**(7/3)*f**(14/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_47(f):
    return 31017210790.9953*(0.024105978576213*pi**(2/3)*f**(2/3) - 0.00315231420649041*pi**2*f + 1)/(pi**2*f**(13/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_55(f):
    return 2047654949.68366*(pi**(2/3)*f**(2/3) - 0.722774823465934*pi**2*f)**2/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_56(f):
    return -41893658675.5658*(0.0147197441049309*pi**(2/3)*f**(2/3) - 0.0106390604469051*pi**2*f)/(pi**(7/3)*f**(14/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_57(f):
    return 8271256210.93209*(0.0147197441049309*pi**(2/3)*f**(2/3) - 0.0106390604469051*pi**2*f)/(pi**2*f**(13/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_66(f):
    return 185712152.913065/(pi**(4/3)*f**(11/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_67(f):
    return -36665997.8333108/(pi*f**(10/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_77(f):
    return 7239135.2748018/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))


Gamma = np.empty((7,7))
f_max = (6**(3/2)*pi*20*(M_sun*G/c**3).value)**(-1)
for i in range(1,8):
    for j in range(i,8):
        exec(f'quad_prod_ij = 2*quad(prod_{i}{j}, 10, f_max)[0]')
        Gamma[i-1,j-1] = quad_prod_ij
        Gamma[j-1,i-1] = quad_prod_ij

with open(__file__[:-3]+'_Gamma.txt', 'w') as f:
    for i in range(7):
        for j in range(7):
            f.write(str(Gamma[i,j])+' ')
        f.write('\n')

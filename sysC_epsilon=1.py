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
    return -1091431799.48701*(4.20800062139669e-5*pi**(4/3)*f**(4/3) + 0.0082441287072856*pi**(2/3)*f**(2/3) - 0.000630462841298082*pi**2*f + 1)/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_25(f):
    return -291048479.863202*(0.000186060467725383*pi**(4/3)*f**(4/3) + 0.00503408167213358*pi**(2/3)*f**(2/3) - 0.00212781208938103*pi**2*f)/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_26(f):
    return 258040.197520371*pi**(1/3)/(f**2*(f**2/2450 + 2 + 24010000/f**4))

def prod_27(f):
    return -29793.4329677736*pi**(2/3)/(f**(5/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_33(f):
    return 51370.1092733744/(f**(7/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_34(f):
    return 38200112982.0452*(4.20800062139669e-5*pi**(4/3)*f**(4/3) + 0.0082441287072856*pi**(2/3)*f**(2/3) - 0.000630462841298082*pi**2*f + 1)/(pi**(5/3)*f**4*(f**2/2450 + 2 + 24010000/f**4))

def prod_35(f):
    return 10186696795.2121*(0.000186060467725383*pi**(4/3)*f**(4/3) + 0.00503408167213358*pi**(2/3)*f**(2/3) - 0.00212781208938103*pi**2*f)/(pi**(5/3)*f**4*(f**2/2450 + 2 + 24010000/f**4))

def prod_36(f):
    return -9031406.91321299/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_37(f):
    return 1042770.15387208/(pi**(1/3)*f**(8/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_44(f):
    return 2.84065705228578e+16*(4.20800062139669e-5*pi**(4/3)*f**(4/3) + 0.0082441287072856*pi**(2/3)*f**(2/3) - 0.000630462841298082*pi**2*f + 1)**2/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_45(f):
    return 7.57508547276208e+15*(0.000186060467725383*pi**(4/3)*f**(4/3) + 0.00503408167213358*pi**(2/3)*f**(2/3) - 0.00212781208938103*pi**2*f)*(4.20800062139669e-5*pi**(4/3)*f**(4/3) + 0.0082441287072856*pi**(2/3)*f**(2/3) - 0.000630462841298082*pi**2*f + 1)/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_46(f):
    return -6715982686265.72*(4.20800062139669e-5*pi**(4/3)*f**(4/3) + 0.0082441287072856*pi**(2/3)*f**(2/3) - 0.000630462841298082*pi**2*f + 1)/(pi**(7/3)*f**(14/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_47(f):
    return 775430269774.884*(4.20800062139669e-5*pi**(4/3)*f**(4/3) + 0.0082441287072856*pi**(2/3)*f**(2/3) - 0.000630462841298082*pi**2*f + 1)/(pi**2*f**(13/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_55(f):
    return 51191373742.0915*(0.0369601607290819*pi**(4/3)*f**(4/3) + pi**(2/3)*f**(2/3) - 0.422681280909613*pi**2*f)**2/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_56(f):
    return -1790928716337.52*(0.000186060467725383*pi**(4/3)*f**(4/3) + 0.00503408167213358*pi**(2/3)*f**(2/3) - 0.00212781208938103*pi**2*f)/(pi**(7/3)*f**(14/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_57(f):
    return 206781405273.302*(0.000186060467725383*pi**(4/3)*f**(4/3) + 0.00503408167213358*pi**(2/3)*f**(2/3) - 0.00212781208938103*pi**2*f)/(pi**2*f**(13/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_66(f):
    return 1587816572.43443/(pi**(4/3)*f**(11/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_67(f):
    return -183329989.166554/(pi*f**(10/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_77(f):
    return 21167359.952843/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))


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

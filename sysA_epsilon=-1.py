from math import pi

import numpy as np
from astropy.constants import G, M_sun, c
from scipy.integrate import quad


def prod_11(f):
    return 50808.2057935719/(f**(7/3)*(f**2/2450 + 2 + 24010000/f**4))

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
    return 41.4760863620995*pi**2/(f**(1/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_23(f):
    return -1451.66302267348*pi/(f**(4/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_24(f):
    return -51821484063.0761*(-1.90127044062539e-6*pi**(4/3)*f**(4/3) + 0.00175238212206869*pi**(2/3)*f**(2/3) - 6.1785358447212e-5*pi**2*f + 1)/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_25(f):
    return -13819062416.8203*(-8.40663534260103e-6*pi**(4/3)*f**(4/3) + 0.00107005058223856*pi**(2/3)*f**(2/3) - 0.000208525584759341*pi**2*f)/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_26(f):
    return 1200678.36303886*pi**(1/3)/(f**2*(f**2/2450 + 2 + 24010000/f**4))

def prod_27(f):
    return -63914.867360584*pi**(2/3)/(f**(5/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_33(f):
    return 50808.2057935719/(f**(7/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_34(f):
    return 1813751942207.66*(-1.90127044062539e-6*pi**(4/3)*f**(4/3) + 0.00175238212206869*pi**(2/3)*f**(2/3) - 6.1785358447212e-5*pi**2*f + 1)/(pi**(5/3)*f**4*(f**2/2450 + 2 + 24010000/f**4))

def prod_35(f):
    return 483667184588.71*(-8.40663534260103e-6*pi**(4/3)*f**(4/3) + 0.00107005058223856*pi**(2/3)*f**(2/3) - 0.000208525584759341*pi**2*f)/(pi**(5/3)*f**4*(f**2/2450 + 2 + 24010000/f**4))

def prod_36(f):
    return -42023742.7063602/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))

def prod_37(f):
    return 2237020.35762044/(pi**(1/3)*f**(8/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_44(f):
    return 6.47473386725707e+19*(-1.90127044062539e-6*pi**(4/3)*f**(4/3) + 0.00175238212206869*pi**(2/3)*f**(2/3) - 6.1785358447212e-5*pi**2*f + 1)**2/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_45(f):
    return 1.72659569793522e+19*(-8.40663534260103e-6*pi**(4/3)*f**(4/3) + 0.00107005058223856*pi**(2/3)*f**(2/3) - 0.000208525584759341*pi**2*f)*(-1.90127044062539e-6*pi**(4/3)*f**(4/3) + 0.00175238212206869*pi**(2/3)*f**(2/3) - 6.1785358447212e-5*pi**2*f + 1)/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_46(f):
    return -1.50016407314543e+15*(-1.90127044062539e-6*pi**(4/3)*f**(4/3) + 0.00175238212206869*pi**(2/3)*f**(2/3) - 6.1785358447212e-5*pi**2*f + 1)/(pi**(7/3)*f**(14/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_47(f):
    return 79857179662610.4*(-1.90127044062539e-6*pi**(4/3)*f**(4/3) + 0.00175238212206869*pi**(2/3)*f**(2/3) - 6.1785358447212e-5*pi**2*f + 1)/(pi**2*f**(13/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_55(f):
    return 5271910176120.45*(-0.00785629715280774*pi**(4/3)*f**(4/3) + pi**(2/3)*f**(2/3) - 0.194874511747943*pi**2*f)**2/(pi**(10/3)*f**(17/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_56(f):
    return -400043752838781.0*(-8.40663534260103e-6*pi**(4/3)*f**(4/3) + 0.00107005058223856*pi**(2/3)*f**(2/3) - 0.000208525584759341*pi**2*f)/(pi**(7/3)*f**(14/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_57(f):
    return 21295247910029.4*(-8.40663534260103e-6*pi**(4/3)*f**(4/3) + 0.00107005058223856*pi**(2/3)*f**(2/3) - 0.000208525584759341*pi**2*f)/(pi**2*f**(13/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_66(f):
    return 34758065620.8449/(pi**(4/3)*f**(11/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_67(f):
    return -1850251676.26417/(pi*f**(10/3)*(f**2/2450 + 2 + 24010000/f**4))

def prod_77(f):
    return 98493146.9680318/(pi**(2/3)*f**3*(f**2/2450 + 2 + 24010000/f**4))


Gamma = np.empty((7,7))
f_max = (6**(3/2)*pi*2.8*(M_sun*G/c**3).value)**(-1)
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

import math as m

from astropy.constants import G, M_sun, c
from scipy.integrate import quad


def I_7(f):
    return f**(-7/3)/(f**(-4)+2+2*f**2)


print(quad(I_7, 1/7, (6**(3/2)*m.pi*2.8*(M_sun*G/c**3).value*70)**(-1))[0])
# print(quad(I_7, 1/7, (6**(3/2)*m.pi*11.4*(M_sun*G/c**3).value*70)**(-1))[0])
# print(quad(I_7, 1/7, (6**(3/2)*m.pi*20*(M_sun*G/c**3).value*70)**(-1))[0])

import math as m

from astropy.constants import G, M_sun, c
from sympy import *


f = symbols('f', positive=True)
PARAMS = ('A_4', 'B_4', 'C_4', 'A_5', 'B_5', 'C_5', 'A', 'M', 'S_0', 'SNR',
          'beta', 'epsilon', 'eta', 'sigma', 'f_0', 'h_F', 'v')
for param in PARAMS:
    exec(f'{param}=symbols("{param}")')
prod = symbols("prod")
FRAC = Rational


S = FRAC(1,5)*S_0*((f_0/f)**4+2+2*(f/f_0)**2)

h_F_1 = (1)*h_F
h_F_2 = (2*pi*I*(f/f_0))*h_F
h_F_3 = (-I)*h_F
h_F_4 = (FRAC(5,128)*I*(pi*M*f)**(-FRAC(5,3))
         *(1+A_4*v**2-B_4*v**3+C_4*v**4))*h_F
h_F_5 = (FRAC(1,96)*I*(pi*M*f)**(-FRAC(5,3))
         *(A_5*v**2-B_5*v**3+C_5*v**4))*h_F
h_F_6 = (FRAC(3,32)*I*(pi*M*f)**(-FRAC(2,3))
         *eta**(-FRAC(3,5)))*h_F
h_F_7 = (FRAC(15,64)*I*(pi*M*f)**(-FRAC(1,3))
         *eta**(-FRAC(4,5)))*h_F

h_F_4 = h_F_4.subs(A_4, FRAC(4,3)*(FRAC(743,336)+FRAC(11,4)*eta))
h_F_4 = h_F_4.subs(B_4, FRAC(8,5)*(4*pi-beta))
h_F_4 = h_F_4.subs(C_4, 2*epsilon*(FRAC(3058673,1016064)+FRAC(5429,1008)*eta
                        +FRAC(617,144)*eta**2-sigma))
h_F_5 = h_F_5.subs(A_5, FRAC(743,168)-FRAC(33,4)*eta)
h_F_5 = h_F_5.subs(B_5, FRAC(27,5)*(4*pi-beta))
h_F_5 = h_F_5.subs(C_5, 18*epsilon*(FRAC(3058673,1016064)-FRAC(5429,4032)*eta
                        -FRAC(617,96)*eta**2-sigma))
for n in range(1,8):
    exec(f'h_F_{n}=h_F_{n}.subs(v,(pi*M/eta**FRAC(3,5)*f)**FRAC(1,3))')

S = S.subs(f_0, 10)
for n in range(1,8):
    exec(f'h_F_{n}=h_F_{n}.subs(f_0,70)')
    exec(f'h_F_{n}=h_F_{n}.subs(M,0.25**(3/5)*2.8*(M_sun*G/c**3).value)')
    exec(f'h_F_{n}=h_F_{n}.subs(eta,0.25)')
    exec(f'h_F_{n}=h_F_{n}.subs(beta,0)')
    exec(f'h_F_{n}=h_F_{n}.subs(sigma,0)')
    exec(f'h_F_{n}=h_F_{n}.subs(epsilon,1)')

SNR = (20*A**2*S_0**(-1)*f_0**(-FRAC(4,3))*0.28390566697806996)
print(SNR)

import math as m

from astropy.constants import G, M_sun, c
from sympy import *


f = symbols('f', positive=True)
PARAMS = ('A_4', 'B_4', 'C_4', 'A_5', 'B_5', 'C_5', 'A', 'M', 'S_0', 'SNR',
          'beta', 'epsilon', 'eta', 'sigma', 'f_0', 'h_F', 'v')
for param in PARAMS:
    exec(f'{param}=symbols("{param}")')
prod = symbols("prod")

conj = conjugate
frac = Rational


S = frac(1,5)*S_0*((f_0/f)**4+2+2*(f/f_0)**2)
# SNR = (20*A**2*S_0**(-1)*f_0**(-frac(4,3))*I_7)

h_F_1 = (1)*h_F
h_F_2 = (2*pi*I*(f/f_0))*h_F
h_F_3 = (-I)*h_F
h_F_4 = (-frac(5,128)*I*(pi*M*f)**(-frac(5,3))
         *(1+A_4*v**2-B_4*v**3+C_4*v**4))*h_F
h_F_5 = (-frac(1,96)*I*(pi*M*f)**(-frac(5,3))
         *(A_5*v**2-B_5*v**3+C_5*v**4))*h_F
h_F_6 = (frac(3,32)*I*(pi*M*f)**(-frac(2,3))
         *eta**(-frac(3,5)))*h_F
h_F_7 = (-frac(15,64)*I*(pi*M*f)**(-frac(1,3))
         *eta**(-frac(4,5)))*h_F

h_F_4 = h_F_4.subs(A_4, frac(4,3)*(frac(743,336)+frac(11,4)*eta))
h_F_4 = h_F_4.subs(B_4, frac(8,5)*(4*pi-beta))
h_F_4 = h_F_4.subs(C_4, 2*epsilon*(frac(3058673,1016064)+frac(5429,1008)*eta
                        +frac(617,144)*eta**2-sigma))
h_F_5 = h_F_5.subs(A_5, frac(743,168)-frac(33,4)*eta)
h_F_5 = h_F_5.subs(B_5, frac(27,5)*(4*pi-beta))
h_F_5 = h_F_5.subs(C_5, 18*epsilon*(frac(3058673,1016064)-frac(5429,4032)*eta
                        -frac(617,96)*eta**2-sigma))
for n in range(1,8):
    exec(f'h_F_{n}=h_F_{n}.subs(v,(pi*M/eta**frac(3,5)*f)**frac(1,3))')

eta_value = (1.4*1.4)/(1.4+1.4)**2
M_value = (1.4*1.4)*(M_sun*G/c**3).value
# eta_value = (1.4*10)/(1.4+10)**2
# M_value = (1.4*10)*(M_sun*G/c**3).value
# eta_value = (10*10)/(10+10)**2
# M_value = (10*10)*(M_sun*G/c**3).value

epsilon_value = 1
# epsilon_value = 0
# epsilon_value = -1

S = S.subs(f_0, 70)
# SNR = SNR.subs(f_0, 70)
for n in range(1,8):
    exec(f'h_F_{n}=h_F_{n}.subs(f_0,70)')
    exec(f'h_F_{n}=h_F_{n}.subs(M,eta_value**(3/5)*M_value)')    
    exec(f'h_F_{n}=h_F_{n}.subs(eta,eta_value)')
    exec(f'h_F_{n}=h_F_{n}.subs(beta,0)')
    exec(f'h_F_{n}=h_F_{n}.subs(sigma,0)')
    exec(f'h_F_{n}=h_F_{n}.subs(epsilon,epsilon_value)')

I_7 = 0.28390096292585254
# I_7 = 0.2834082747325775
# I_7 = 0.2807955590004223
for i in range(1,8):
    for j in range(i,8):
        exec(f'prod = (conj(h_F_{i})*h_F_{j}+h_F_{i}*conj(h_F_{j}))/S')
        prod = prod.subs(h_F*conj(h_F), (A*f**(-frac(7,6)))**2)
        prod = prod.subs(A**2*S_0**(-1), SNR/(20*70**(-4/3)*I_7))
        prod = prod.subs(SNR, 10**2)
        print(f'i={i}, j={j}')
        print(prod)
        print('\n')

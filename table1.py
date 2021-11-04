from astropy.constants import G, M_sun, c
from sympy import *


for symbol in ('N_newtonian', 'N_1pn', 'N_tail',
               'N_spin_orbit', 'N_2pn', 'N_spin_spin',
               'F', 'eta', 'M', 'F_min', 'F_max', 'beta', 'sigma'):
    exec(f'{symbol}=symbols("{symbol}")')
frac = Rational


# N=\int_{f_\text{min}}^{f_\text{max}}\frac{F}{\dot{F}}dF
N_newtonian = integrate(
    F
    *(5*pi*(eta**frac(3,5)*M)**2/96)
    *(pi*eta**frac(3,5)*M*F)**frac(-11,3),
    (F, F_min, F_max)
)
N_1pn = integrate(
    F
    *(5*pi*(eta**frac(3,5)*M)**2/96)
    *(pi*eta**frac(3,5)*M*F)**frac(-11,3)
    *(frac(743,336)+frac(11,4)*eta)
    *(pi*M*F)**frac(2,3),
    (F, F_min, F_max)
)
N_tail = integrate(
    F
    *(5*pi*(eta**frac(3,5)*M)**2/96)
    *(pi*eta**frac(3,5)*M*F)**frac(-11,3)
    *(-4*pi)
    *(pi*M*F),
    (F, F_min, F_max)
)
N_spin_orbit = integrate(
    F
    *(5*pi*(eta**frac(3,5)*M)**2/96)
    *(pi*eta**frac(3,5)*M*F)**frac(-11,3)
    *(beta)
    *(pi*M*F),
    (F, F_min, F_max)
)
N_2pn = integrate(
    F
    *(5*pi*(eta**frac(3,5)*M)**2/96)
    *(pi*eta**frac(3,5)*M*F)**frac(-11,3)
    *(frac(3058673,1016064)+frac(5429,1008)*eta+frac(617,144)*eta**2) # FIXME: please check this expression
    *(pi*M*F)**frac(4,3),
    (F, F_min, F_max)
)
N_spin_spin = integrate(
    F
    *(5*pi*(eta**frac(3,5)*M)**2/96)
    *(pi*eta**frac(3,5)*M*F)**frac(-11,3)
    *(-sigma)
    *(pi*M*F)**frac(4,3),
    (F, F_min, F_max)
)

for N_ in (N_newtonian, N_1pn, N_tail,
           N_spin_orbit, N_2pn, N_spin_spin):
    N = N_
    N = N.subs(eta, (1.4*1.4)/(1.4+1.4)**2)
    N = N.subs(M, (1.4+1.4)*(M_sun*G/c**3).value)
    N = N.subs(F_min, 10)
    N = N.subs(F_max, 1000)
    print(str(N.evalf())+' ', end='')
    N = N_
    N = N.subs(eta, (1.4*10)/(1.4+10)**2)
    N = N.subs(M, (1.4+10)*(M_sun*G/c**3).value)
    N = N.subs(F_min, 10)
    N = N.subs(F_max, 360)
    print(str(N.evalf())+' ', end='')
    N = N_
    N = N.subs(eta, (10*10)/(10+10)**2)
    N = N.subs(M, (10+10)*(M_sun*G/c**3).value)
    N = N.subs(F_min, 10)
    N = N.subs(F_max, 190)
    print(str(N.evalf())+' ', end='')
    print('\n')

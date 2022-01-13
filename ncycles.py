from math import pi

from astropy.constants import G, M_sun, c
from scipy.integrate import quad 


M_SUN = (M_sun*G/c**3).value
TERMS = ('newtonian', '1pn', 'tail', 'spin orbit', '2pn', 'spin spin')


def cycles_of_binary_GW(m1, m2, f_min, f_max, terms=TERMS):
    '''
    The accumulated number of cycles of gravitational waves from binaries.
    '''

    m1, m2 = m1*M_SUN, m2*M_SUN

    M = m1+m2
    mu = m1*m2/(m1+m2)
    eta = mu/M
    calM = eta**(3/5)*M

    params_of_terms = {
        'newtonian': (1, 0),
        '1pn': ((743/336)+(11/4)*eta, 2/3),
        'tail': (-4*pi, 1),
        'spin orbit': (1, 1), # set beta = 1
        '2pn': ((3058673/1016064)+(5429/1008)*eta+(617/144)*eta**2, 4/3),
        'spin spin': (-1, 4/3), # set sigma = 1
    }

    cycles_of_terms = {}
    for term in terms:
        def cycles_per_unit_freq(f):
            params = params_of_terms[term]
            N_f = (
               f*(5*pi*calM**2/96)*(pi*calM*f)**(-11/3)
               *params[0]*(pi*M*f)**params[1] 
            )
            return N_f
        cycles_of_terms[term] = quad(cycles_per_unit_freq, f_min, f_max)[0]
    
    return cycles_of_terms


if (__name__ == '__main__'):
    print(cycles_of_binary_GW(1.4, 1.4, 10, 1000))
    print(cycles_of_binary_GW(1.4, 10, 10, 360))
    print(cycles_of_binary_GW(10, 10, 10, 190))

import numpy as np

Msun = 1.989e30 # mass of sun, kg
G = 6.67408e-11 * Msun # m3 Msun-1 s-2
c = 299792458 # m/s
M2T = G/c**3

# Poisson & Will 1995

def cycleGWfromdFdtPN(m1, m2, f_min, f_max=None, z=0):
    m1, m2 = m1 * (1+z), m2 * (1+z)
    M = m1 + m2
    eta = m1*m2 / M**2
    chirpM = eta ** (3/5) * M
    
    dfdtPNList = {
        'Newtonian':[0, -1], # When we use delta N formula, should use -1 (not 1) 
        '1PN':[1, -(743/336 + 11/4 * eta)], 
        'Tail': [1.5, 4 * np.pi], # set beta = 0, also Tail term
        'SO with beta=1': [1.5, -1], # set beta = 1
        '2PN': [2, (34103/18144 + 13661/2016 * eta + 59/18 * eta**2)], # set sigma=0
        'SS with sigma=1': [2, 1] # set sigma=1
        }

    if f_max == None:
        f_max = 1/np.pi * M * M2T * 6**(2/3)
        
    f = np.logspace(np.log10(f_min), np.log10(f_max), 1000, endpoint=True)
    u = np.pi * M * M2T * f

    dfdtNewtonian = 96/(5*np.pi*chirpM**2 * M2T**2) * (np.pi * chirpM * M2T * f)**(11/3)
    
    results = {}
    for mode in dfdtPNList.keys():
        pn, term = dfdtPNList[mode]
        DeltaN = np.trapz(-f * term * u**(2/3 * pn)/dfdtNewtonian, f)
        results[mode] = DeltaN

    return results

if __name__ == '__main__':
    binaries = [(1.4, 1.4, 10, 1000),
           (1.4, 10, 10, 360),
           (10, 10, 10, 190)]
    for binary in binaries:
        print('\n m1=%.2f, m2=%.2f, f_min=%.2f, f_max=%.2f\n'%binary)
        res = cycleGWfromdFdtPN(*binary)
        for mode in res.keys():
            print('%s: %.2f'%(mode, res[mode]))
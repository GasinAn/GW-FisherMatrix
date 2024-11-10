from numpy import *
from scipy.integrate import *


G = 6.6743e-11
c = 299792458
z = 1
beta = 0
phi_c = 0
t_c = 0
M_sun = 1.988409870698051e+30*(G/c**3)
R = 149597870700/c
T = sqrt(R**3/(M_sun/(4*pi**2)))
ln_D_L = log(z/(75*1e3/(1e6*(R*c)/(pi/180/60/60))))


def input2param(M_1, M_2):
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    return log(cal_M), log(mu)

def h(Alpha, f, ln_cal_M, ln_mu, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    alpha_0 = 0
    bar_phi_0 = 0

    cal_M, mu = exp(ln_cal_M), exp(ln_mu)
    M = cal_M**(5/2)/mu**(3/2)
    eta = mu/M

    x = (pi*M*(1+z)*f)**(2/3)
    t = t_c-5*(8*pi*f)**(-8/3)*(cal_M*(1+z))**(-5/3)*(
        1+(4/3)*((743/336)+(11/4)*eta)*x-(8/5)*(4*pi-beta)*x**(3/2)
    )
    Psi = 2*pi*f*t_c-phi_c-pi/4+(3/4)*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(
        1+(20/9)*((743/336)+(11/4)*eta)*x-4*(4*pi-beta)*x**(3/2)
    )

    bar_theta_S, bar_theta_L = arccos(bar_mu_S), arccos(bar_mu_L)
    bar_phi = bar_phi_0+2*pi*t/T

    theta_S = arccos(
              ((sqrt(1)/2)*cos(bar_theta_S)
              -(sqrt(3)/2)*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)))

    phi_S = alpha_0+2*pi*t/T+arctan(
            (sin(bar_theta_S)*cos(bar_phi-bar_phi_S)+sqrt(3)*cos(bar_theta_S))
            /(2*sin(bar_theta_S)*sin(bar_phi-bar_phi_S)))

    Lz = ((sqrt(1)/2)*cos(bar_theta_L)
         -(sqrt(3)/2)*sin(bar_theta_L)*cos(bar_phi-bar_phi_L))
    Ln = (cos(bar_theta_L)*cos(bar_theta_S)
         +sin(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_L-bar_phi_S))
    zn = cos(theta_S)
    nLz = ((1/2)*sin(bar_theta_L)*sin(bar_theta_S)*sin(bar_phi_L-bar_phi_S)
            -(sqrt(3)/2)*cos(bar_phi)
                        *(cos(bar_theta_L)*sin(bar_theta_S)*sin(bar_phi_S)
                         -cos(bar_theta_S)*sin(bar_theta_L)*sin(bar_phi_L))
            -(sqrt(3)/2)*sin(bar_phi)
                        *(cos(bar_theta_S)*sin(bar_theta_L)*cos(bar_phi_L)
                         -cos(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_S)))
    psi_S = arctan((Lz-Ln*zn)/nLz)

    F_p = (cos(2*(phi_S-(Alpha-1)*(pi/4)))*cos(2*psi_S)*(1+cos(theta_S)**2)/2
          -sin(2*(phi_S-(Alpha-1)*(pi/4)))*sin(2*psi_S)*cos(theta_S))
    F_t = (cos(2*(phi_S-(Alpha-1)*(pi/4)))*sin(2*psi_S)*(1+cos(theta_S)**2)/2
          +sin(2*(phi_S-(Alpha-1)*(pi/4)))*cos(2*psi_S)*cos(theta_S))
    Lambda_p, Lambda_t = 1+Ln**2, -2*Ln
    Lambda = sqrt((Lambda_p*F_p)**2+(Lambda_t*F_t)**2)

    cal_A = (5/96)**(1/2)*pi**(-2/3)*(cal_M*(1+z))**(5/6)*(exp(ln_D_L))**(-1)
    A = (sqrt(3)/2)*Lambda*cal_A*f**(-7/6)

    varphi_p = arctan(-(Lambda_t*F_t)/(Lambda_p*F_p))
    varphi_D = 2*pi*f*R*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)

    return A*exp(1j*(Psi-varphi_p-varphi_D))

def h_0(**kwargs):
    ln_cal_M_p = kwargs['ln_cal_M']+kwargs['ln_cal_M']*1e-5
    ln_cal_M_n = kwargs['ln_cal_M']-kwargs['ln_cal_M']*1e-5
    kwargs_p, kwargs_n = kwargs.copy(), kwargs.copy()
    kwargs_p['ln_cal_M'], kwargs_n['ln_cal_M'] = ln_cal_M_p, ln_cal_M_n
    return (h(**kwargs_p)-h(**kwargs_n))/(kwargs['ln_cal_M']*2e-5)

def h_1(**kwargs):
    ln_mu_p = kwargs['ln_mu']+kwargs['ln_mu']*1e-5
    ln_mu_n = kwargs['ln_mu']-kwargs['ln_mu']*1e-5
    kwargs_p, kwargs_n = kwargs.copy(), kwargs.copy()
    kwargs_p['ln_mu'], kwargs_n['ln_mu'] = ln_mu_p, ln_mu_n
    return (h(**kwargs_p)-h(**kwargs_n))/(kwargs['ln_mu']*2e-5)

def h_2(**kwargs):
    f, ln_cal_M, ln_mu = kwargs['f'], kwargs['ln_cal_M'], kwargs['ln_mu']

    cal_M, mu = exp(ln_cal_M), exp(ln_mu)
    M = cal_M**(5/2)/mu**(3/2)

    return (3*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)
           *(pi*M*(1+z)*f)*h(**kwargs))

def h_3(**kwargs):
    return -1j*h(**kwargs)

def h_4(**kwargs):
    return 2*pi*1j*kwargs['f']*h(**kwargs)

def h_5(**kwargs):
    return -h(**kwargs)

def h_6(**kwargs):
    bar_mu_S_p = kwargs['bar_mu_S']+1e-5
    bar_mu_S_n = kwargs['bar_mu_S']-1e-5
    kwargs_p, kwargs_n = kwargs.copy(), kwargs.copy()
    kwargs_p['bar_mu_S'], kwargs_n['bar_mu_S'] = bar_mu_S_p, bar_mu_S_n
    return (h(**kwargs_p)-h(**kwargs_n))/2e-5

def h_7(**kwargs):
    bar_phi_S_p = kwargs['bar_phi_S']+1e-5
    bar_phi_S_n = kwargs['bar_phi_S']-1e-5
    kwargs_p, kwargs_n = kwargs.copy(), kwargs.copy()
    kwargs_p['bar_phi_S'], kwargs_n['bar_phi_S'] = bar_phi_S_p, bar_phi_S_n
    return (h(**kwargs_p)-h(**kwargs_n))/2e-5

def h_8(**kwargs):
    bar_mu_L_p = kwargs['bar_mu_L']+1e-5
    bar_mu_L_n = kwargs['bar_mu_L']-1e-5
    kwargs_p, kwargs_n = kwargs.copy(), kwargs.copy()
    kwargs_p['bar_mu_L'], kwargs_n['bar_mu_L'] = bar_mu_L_p, bar_mu_L_n
    return (h(**kwargs_p)-h(**kwargs_n))/2e-5

def h_9(**kwargs):
    bar_phi_L_p = kwargs['bar_phi_L']+1e-5
    bar_phi_L_n = kwargs['bar_phi_L']-1e-5
    kwargs_p, kwargs_n = kwargs.copy(), kwargs.copy()
    kwargs_p['bar_phi_L'], kwargs_n['bar_phi_L'] = bar_phi_L_p, bar_phi_L_n
    return (h(**kwargs_p)-h(**kwargs_n))/2e-5

def S_n(f):
    f = abs(f)

    alpha = 10**(-22.79)*(f/1e-3)**(-7/3)
    beta  = 10**(-24.54)*(f/1e-3)
    gamma = 10**(-23.04)
    S_n_in = 5.049e+5*(alpha**2+beta**2+gamma**2)

    if (f <= 10**(-3.15)):
        S_n_co = 10**(-42.685)*f**(-1.9)
    elif (f <= 10**(-2.75)):
        S_n_co = 10**(-60.325)*f**(-7.5)
    else:
        S_n_co = 10**(-46.850)*f**(-2.6)

    return S_n_in+S_n_co

def f_max(M_1, M_2):
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    return (3**(3/2)*pi*M*(1+z))**(-1)

def signal2noise_I(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    ln_cal_M, ln_mu = input2param(M_1, M_2)
    kwargs = {'ln_cal_M': ln_cal_M, 'ln_mu': ln_mu,
              'bar_mu_S': bar_mu_S, 'bar_phi_S': bar_phi_S,
              'bar_mu_L': bar_mu_L, 'bar_phi_L': bar_phi_L}
    def func_integrated(f):
        kwargs['Alpha'], kwargs['f'] = 1, f
        h_I = h(**kwargs)
        return (abs(h_I)**2)/S_n(f)
    square_rho = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return sqrt(square_rho)

def signal2noise(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    ln_cal_M, ln_mu = input2param(M_1, M_2)
    kwargs = {'ln_cal_M': ln_cal_M, 'ln_mu': ln_mu,
              'bar_mu_S': bar_mu_S, 'bar_phi_S': bar_phi_S,
              'bar_mu_L': bar_mu_L, 'bar_phi_L': bar_phi_L}
    def func_integrated(f):
        kwargs['Alpha'], kwargs['f'] = 1, f
        h_I = h(**kwargs)
        kwargs['Alpha'], kwargs['f'] = 2, f
        h_II = h(**kwargs)
        return (abs(h_I)**2+abs(h_II)**2)/S_n(f)
    square_rho = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return sqrt(square_rho)

def Fisher_matrix_I(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    ln_cal_M, ln_mu = input2param(M_1, M_2)
    kwargs = {'ln_cal_M': ln_cal_M, 'ln_mu': ln_mu,
              'bar_mu_S': bar_mu_S, 'bar_phi_S': bar_phi_S,
              'bar_mu_L': bar_mu_L, 'bar_phi_L': bar_phi_L}
    Gamma = empty((10,10))
    for i in range(10):
        for j in range(10):
            def func_integrated(f):
                kwargs['Alpha'], kwargs['f'] = 1, f
                h_I_i = eval(f'h_{i}(**kwargs)')
                h_I_j = eval(f'h_{j}(**kwargs)')
                return (h_I_i.real*h_I_j.real+h_I_i.imag*h_I_j.imag)/S_n(f)
            Gamma[i,j] = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return Gamma

def Fisher_matrix(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    ln_cal_M, ln_mu = input2param(M_1, M_2)
    kwargs = {'ln_cal_M': ln_cal_M, 'ln_mu': ln_mu,
              'bar_mu_S': bar_mu_S, 'bar_phi_S': bar_phi_S,
              'bar_mu_L': bar_mu_L, 'bar_phi_L': bar_phi_L}
    Gamma = empty((10,10))
    for i in range(10):
        for j in range(10):
            def func_integrated(f):
                kwargs['Alpha'], kwargs['f'] = 1, f
                h_I_i = eval(f'h_{i}(**kwargs)')
                h_I_j = eval(f'h_{j}(**kwargs)')
                kwargs['Alpha'], kwargs['f'] = 2, f
                h_II_i = eval(f'h_{i}(**kwargs)')
                h_II_j = eval(f'h_{j}(**kwargs)')
                return (h_I_i.real*h_I_j.real+h_I_i.imag*h_I_j.imag
                       +h_II_i.real*h_II_j.real+h_II_i.imag*h_II_j.imag)/S_n(f)
            Gamma[i,j] = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return Gamma

def Culter(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = (M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    rho_I, rho = signal2noise_I(*args), signal2noise(*args)
    Gamma_I, Gamma = Fisher_matrix_I(*args), Fisher_matrix(*args)
    Sigma_I, Sigma = linalg.inv(Gamma_I), linalg.inv(Gamma)
    Delta_Omega_S_I = 2*pi*(sqrt(Sigma_I[6,6])*sqrt(Sigma_I[7,7])-Sigma_I[6,7])
    Delta_Omega_S = 2*pi*(sqrt(Sigma[6,6])*sqrt(Sigma[7,7])-Sigma[6,7])
    Delta_ln_D_L, Delta_ln_mu = sqrt(Sigma[5,5]), sqrt(Sigma[1,1])
    return rho_I, rho, Delta_Omega_S_I, Delta_Omega_S, Delta_ln_D_L, Delta_ln_mu


print(Culter(1e7,1e7,0.3,5.0,0.8,2.0))
print(Culter(1e6,1e6,-0.8,1.0,0.5,3.0))
print(Culter(1e5,1e5,0.9,2.0,-0.8,5.0))
print(Culter(1e4,1e4,-0.1,3.0,-0.9,6.0))

#print(signal2noise_I(1e7,1e7,0.3,5.0,0.8,2.0))
#print(signal2noise(1e7,1e7,0.3,5.0,0.8,2.0))
#print(signal2noise_I(1e6,1e6,-0.8,1.0,0.5,3.0))
#print(signal2noise(1e6,1e6,-0.8,1.0,0.5,3.0))
#print(signal2noise_I(1e5,1e5,0.9,2.0,-0.8,5.0))
#print(signal2noise(1e5,1e5,0.9,2.0,-0.8,5.0))
#print(signal2noise_I(1e4,1e4,-0.1,3.0,-0.9,6.0))
#print(signal2noise(1e4,1e4,-0.1,3.0,-0.9,6.0))

#print('1e7')
#Gamma=Fisher_matrix(1e7,1e7,0.3,5.0,0.8,2.0)
#print(Gamma)
#print(linalg.inv(Gamma))
#savetxt('1e7_Gamma_n.txt',Gamma)
#savetxt('1e7_Sigma_n.txt',linalg.inv(Gamma))
#print('1e6')
#Gamma=Fisher_matrix(1e6,1e6,-0.8,1.0,0.5,3.0)
#print(Gamma)
#print(linalg.inv(Gamma))
#savetxt('1e6_Gamma_n.txt',Gamma)
#savetxt('1e6_Sigma_n.txt',linalg.inv(Gamma))
#print('1e5')
#Gamma=Fisher_matrix(1e5,1e5,0.9,2.0,-0.8,5.0)
#print(Gamma)
#print(linalg.inv(Gamma))
#savetxt('1e5_Gamma_n.txt',Gamma)
#savetxt('1e5_Sigma_n.txt',linalg.inv(Gamma))
#print('1e4')
#Gamma=Fisher_matrix(1e4,1e4,-0.1,3.0,-0.9,6.0)
#print(Gamma)
#print(linalg.inv(Gamma))
#savetxt('1e4_Gamma_n.txt',Gamma)
#savetxt('1e4_Sigma_n.txt',linalg.inv(Gamma))

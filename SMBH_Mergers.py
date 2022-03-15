from numpy import *
from scipy.integrate import *
from scipy.misc import derivative


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


def func_Lambda_alpha(Lambda_p, Lambda_t, F_p_alpha, F_t_alpha):
    return sqrt((Lambda_p*F_p_alpha)**2+(Lambda_t*F_t_alpha)**2)

def func_cal_A(cal_M):
    return (5/96)**(1/2)*pi**(-2/3)*(cal_M*(1+z))**(5/6)*(exp(ln_D_L))**(-1)

def func_Psi(f, cal_M, mu):
    M = cal_M**(5/2)/mu**(3/2)
    eta = mu/M
    x = (pi*M*(1+z)*f)**(2/3)
    return 2*pi*f*t_c-phi_c-pi/4+(3/4)*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(
        1+(20/9)*((743/336)+(11/4)*eta)*x-4*(4*pi-beta)*x**(3/2)
    )

def func_t(f, cal_M, mu):
    M = cal_M**(5/2)/mu**(3/2)
    eta = mu/M
    x = (pi*M*(1+z)*f)**(2/3)
    return t_c-5*(8*pi*f)**(-8/3)*(cal_M*(1+z))**(-5/3)*(
        1+(4/3)*((743/336)+(11/4)*eta)*x-(8/5)*(4*pi-beta)*x**(3/2)
    )

def func_varphi_p_alpha(Lambda_p, Lambda_t, F_p_alpha, F_t_alpha):
    return arctan(-(Lambda_t*F_t_alpha)/(Lambda_p*F_p_alpha))

def func_varphi_D(f, bar_theta_S, bar_phi, bar_phi_S):
    return 2*pi*f*R*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)

def func_Lz(bar_theta_L, bar_phi, bar_phi_L):
    return ((sqrt(1)/2)*cos(bar_theta_L)
           -(sqrt(3)/2)*sin(bar_theta_L)*cos(bar_phi-bar_phi_L))

def func_Ln(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S):
    return (cos(bar_theta_L)*cos(bar_theta_S)
           +sin(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_L-bar_phi_S))

def func_zn(theta_S):
    return cos(theta_S)

def func_nLz(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S, bar_phi):
    return ((1/2)*sin(bar_theta_L)*sin(bar_theta_S)*sin(bar_phi_L-bar_phi_S)
            -(sqrt(3)/2)*cos(bar_phi)
                        *(cos(bar_theta_L)*sin(bar_theta_S)*sin(bar_phi_S)
                         -cos(bar_theta_S)*sin(bar_theta_L)*sin(bar_phi_L))
            -(sqrt(3)/2)*sin(bar_phi)
                        *(cos(bar_theta_S)*sin(bar_theta_L)*cos(bar_phi_L)
                         -cos(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_S)))

def func_Lambda_p(Ln):
    return 1+Ln**2

def func_Lambda_t(Ln):
    return -2*Ln

def func_F_p_I(theta_S, phi_S, psi_S):
    return (cos(2*phi_S)*cos(2*psi_S)*(1+cos(theta_S)**2)/2
           -sin(2*phi_S)*sin(2*psi_S)*cos(theta_S))

def func_F_t_I(theta_S, phi_S, psi_S):
    return (cos(2*phi_S)*sin(2*psi_S)*(1+cos(theta_S)**2)/2
           +sin(2*phi_S)*cos(2*psi_S)*cos(theta_S))

def func_F_p_II(theta_S, phi_S, psi_S):
    return (sin(2*phi_S)*cos(2*psi_S)*(1+cos(theta_S)**2)/2
           +cos(2*phi_S)*sin(2*psi_S)*cos(theta_S))

def func_F_t_II(theta_S, phi_S, psi_S):
    return (sin(2*phi_S)*sin(2*psi_S)*(1+cos(theta_S)**2)/2
           -cos(2*phi_S)*cos(2*psi_S)*cos(theta_S))

def func_theta_S(bar_theta_S, bar_phi, bar_phi_S):
    return arccos(
           ((sqrt(1)/2)*cos(bar_theta_S)
           -(sqrt(3)/2)*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)))

def func_phi_S(t, bar_theta_S, bar_phi, bar_phi_S):
    alpha_0 = 0
    return alpha_0+2*pi*t/T+arctan(
        (sin(bar_theta_S)*cos(bar_phi-bar_phi_S)+sqrt(3)*cos(bar_theta_S))
        /(2*sin(bar_theta_S)*sin(bar_phi-bar_phi_S)))

def func_psi_S(Lz, Ln, zn, nLz):
    return arctan((Lz-Ln*zn)/nLz)

def func_bar_phi(t):
    bar_phi_0 = 0
    return bar_phi_0+2*pi*t/T

def S_n(f):
    f = abs(f)
    alpha = 10**(-22.79)*(f/1e-3)**(-7/3)
    beta  = 10**(-24.54)*(f/1e-3)
    gamma = 10**(-23.04)
    S_n_in = 5.049e+5*(alpha**2+beta**2+gamma**2)
    if f<=10**(-3.15):
        S_n_co = 10**(-42.685)*f**(-1.9)
    elif f<=10**(-2.75):
        S_n_co = 10**(-60.325)*f**(-7.5)
    else:
        S_n_co = 10**(-46.850)*f**(-2.6)
    return S_n_in+S_n_co

def f_max(M_1, M_2):
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    return (3**(3/2)*pi*M*(1+z))**(-1)

def preparation(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    bar_theta_S, bar_theta_L = arccos(bar_mu_S), arccos(bar_mu_L)
    t = func_t(f, cal_M, mu)
    bar_phi = func_bar_phi(t)
    theta_S = func_theta_S(bar_theta_S, bar_phi, bar_phi_S)
    Lz = func_Lz(bar_theta_L, bar_phi, bar_phi_L)
    Ln = func_Ln(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S)
    zn = func_zn(theta_S)
    nLz = func_nLz(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S, bar_phi)
    phi_S = func_phi_S(t, bar_theta_S, bar_phi, bar_phi_S)
    psi_S = func_psi_S(Lz, Ln, zn, nLz)
    Lambda_p = func_Lambda_p(Ln)
    Lambda_t = func_Lambda_t(Ln)
    return (t, cal_M, mu, Lambda_p, Lambda_t, theta_S, phi_S, psi_S)

def func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    t, cal_M, mu, Lambda_p, Lambda_t, theta_S, phi_S, psi_S = args
    bar_theta_S = arccos(bar_mu_S)
    bar_phi = func_bar_phi(t)
    F_p_I = func_F_p_I(theta_S, phi_S, psi_S)
    F_t_I = func_F_t_I(theta_S, phi_S, psi_S)
    Lambda_I = func_Lambda_alpha(Lambda_p, Lambda_t, F_p_I, F_t_I)
    cal_A = func_cal_A(cal_M)
    Psi = func_Psi(f, cal_M, mu)
    varphi_p_I = func_varphi_p_alpha(Lambda_p, Lambda_t, F_p_I, F_t_I)
    varphi_D = func_varphi_D(f, bar_theta_S, bar_phi, bar_phi_S)
    A = (sqrt(3)/2)*Lambda_I*cal_A*f**(-7/6)
    return A*exp(1j*(Psi-varphi_p_I-varphi_D))

def func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    t, cal_M, mu, Lambda_p, Lambda_t, theta_S, phi_S, psi_S = args
    bar_theta_S = arccos(bar_mu_S)
    bar_phi = func_bar_phi(t)
    F_p_II = func_F_p_II(theta_S, phi_S, psi_S)
    F_t_II = func_F_t_II(theta_S, phi_S, psi_S)
    Lambda_II = func_Lambda_alpha(Lambda_p, Lambda_t, F_p_II, F_t_II)
    cal_A = func_cal_A(cal_M)
    Psi = func_Psi(f, cal_M, mu)
    varphi_p_II = func_varphi_p_alpha(Lambda_p, Lambda_t, F_p_II, F_t_II)
    varphi_D = func_varphi_D(f, bar_theta_S, bar_phi, bar_phi_S)
    A = (sqrt(3)/2)*Lambda_II*cal_A*f**(-7/6)
    return A*exp(1j*(Psi-varphi_p_II-varphi_D))

def partial_0_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_I = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    eta = mu/M
    x = (pi*M*(1+z)*f)**(2/3)
    return (-5/4)*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(
        1+((55/6)*eta)*x+(8*pi-2*beta)*x**(3/2)
    )*h_I

def partial_1_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_I = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    eta = mu/M
    x = (pi*M*(1+z)*f)**(2/3)
    return (3/4)*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(
        ((-3715/756)+(55/6)*eta)*x+(24*pi-6*beta)*x**(3/2)
    )*h_I

def partial_2_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_I = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    return 3*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(pi*M*(1+z)*f)*h_I

def partial_3_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_I = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return -1j*h_I

def partial_4_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_I = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return 2*pi*1j*f*h_I

def partial_5_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_I = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return -h_I

def partial_6_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_mu_S_p = bar_mu_S+1e-5
    bar_mu_S_n = bar_mu_S-1e-5
    h_I_p = func_h_I(f, M_1, M_2, bar_mu_S_p, bar_phi_S, bar_mu_L, bar_phi_L)
    h_I_n = func_h_I(f, M_1, M_2, bar_mu_S_n, bar_phi_S, bar_mu_L, bar_phi_L)
    return (h_I_p-h_I_n)/2e-5

def partial_7_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_phi_S_p = bar_phi_S+1e-5
    bar_phi_S_n = bar_phi_S-1e-5
    h_I_p = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S_p, bar_mu_L, bar_phi_L)
    h_I_n = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S_n, bar_mu_L, bar_phi_L)
    return (h_I_p-h_I_n)/2e-5

def partial_8_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_mu_L_p = bar_mu_L+1e-5
    bar_mu_L_n = bar_mu_L-1e-5
    h_I_p = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L_p, bar_phi_L)
    h_I_n = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L_n, bar_phi_L)
    return (h_I_p-h_I_n)/2e-5

def partial_9_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_phi_L_p = bar_phi_L+1e-5
    bar_phi_L_n = bar_phi_L-1e-5
    h_I_p = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_p)
    h_I_n = func_h_I(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_n)
    return (h_I_p-h_I_n)/2e-5

def partial_0_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_II = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    eta = mu/M
    x = (pi*M*(1+z)*f)**(2/3)
    return (-5/4)*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(
        1+((55/6)*eta)*x+(8*pi-2*beta)*x**(3/2)
    )*h_II

def partial_1_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_II = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    eta = mu/M
    x = (pi*M*(1+z)*f)**(2/3)
    return (3/4)*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(
        ((-3715/756)+(55/6)*eta)*x+(24*pi-6*beta)*x**(3/2)
    )*h_II

def partial_2_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_II = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    M_1, M_2 = M_1*M_sun, M_2*M_sun
    M = M_1+M_2
    cal_M, mu = (M_1*M_2)**(3/5)/M**(1/5), (M_1*M_2)/M
    return 3*1j*(8*pi*f)**(-5/3)*(cal_M*(1+z))**(-5/3)*(pi*M*(1+z)*f)*h_II

def partial_3_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_II = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return -1j*h_II

def partial_4_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_II = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return 2*pi*1j*f*h_II

def partial_5_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    h_II = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return -h_II

def partial_6_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_mu_S_p = bar_mu_S+1e-5
    bar_mu_S_n = bar_mu_S-1e-5
    h_II_p = func_h_II(f, M_1, M_2, bar_mu_S_p, bar_phi_S, bar_mu_L, bar_phi_L)
    h_II_n = func_h_II(f, M_1, M_2, bar_mu_S_n, bar_phi_S, bar_mu_L, bar_phi_L)
    return (h_II_p-h_II_n)/2e-5

def partial_7_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_phi_S_p = bar_phi_S+1e-5
    bar_phi_S_n = bar_phi_S-1e-5
    h_II_p = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S_p, bar_mu_L, bar_phi_L)
    h_II_n = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S_n, bar_mu_L, bar_phi_L)
    return (h_II_p-h_II_n)/2e-5

def partial_8_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_mu_L_p = bar_mu_L+1e-5
    bar_mu_L_n = bar_mu_L-1e-5
    h_II_p = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L_p, bar_phi_L)
    h_II_n = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L_n, bar_phi_L)
    return (h_II_p-h_II_n)/2e-5

def partial_9_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_phi_L_p = bar_phi_L+1e-5
    bar_phi_L_n = bar_phi_L-1e-5
    h_II_p = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_p)
    h_II_n = func_h_II(f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_n)
    return (h_II_p-h_II_n)/2e-5

def signal2noise_I(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_integrated(f):
        args = (f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
        h_I = func_h_I(*args)
        return (abs(h_I)**2)/S_n(f)
    square_rho = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return sqrt(square_rho)

def signal2noise(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_integrated(f):
        args = (f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
        h_I = func_h_I(*args)
        h_II = func_h_II(*args)
        return (abs(h_I)**2+abs(h_II)**2)/S_n(f)
    square_rho = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return sqrt(square_rho)

def Fisher_matrix_I(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    Gamma = empty((10,10))
    for i in range(10):
        for j in range(10):
            def func_integrated(f):
                args = (f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
                h_I_i = eval(f'partial_{i}_h_I(*args)')
                h_I_j = eval(f'partial_{j}_h_I(*args)')
                return (h_I_i.real*h_I_j.real+h_I_i.imag*h_I_j.imag)/S_n(f)
            Gamma[i,j] = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return Gamma

def Fisher_matrix(M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    Gamma = empty((10,10))
    for i in range(10):
        for j in range(10):
            def func_integrated(f):
                args = (f, M_1, M_2, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
                h_I_i = eval(f'partial_{i}_h_I(*args)')
                h_I_j = eval(f'partial_{j}_h_I(*args)')
                h_II_i = eval(f'partial_{i}_h_II(*args)')
                h_II_j = eval(f'partial_{j}_h_II(*args)')
                return (h_I_i.real*h_I_j.real+h_I_i.imag*h_I_j.imag
                       +h_II_i.real*h_II_j.real+h_II_i.imag*h_II_j.imag)/S_n(f)
            Gamma[i,j] = 4*quad(func_integrated, 0, f_max(M_1, M_2))[0]
    return Gamma


print(signal2noise_I(1e7,1e7,0.3,5.0,0.8,2.0))
print(signal2noise(1e7,1e7,0.3,5.0,0.8,2.0))
print(signal2noise_I(1e6,1e6,-0.8,1.0,0.5,3.0))
print(signal2noise(1e6,1e6,-0.8,1.0,0.5,3.0))
print(signal2noise_I(1e5,1e5,0.9,2.0,-0.8,5.0))
print(signal2noise(1e5,1e5,0.9,2.0,-0.8,5.0))
print(signal2noise_I(1e4,1e4,-0.1,3.0,-0.9,6.0))
print(signal2noise(1e4,1e4,-0.1,3.0,-0.9,6.0))

print('1e7')
Gamma=Fisher_matrix(1e7,1e7,0.3,5.0,0.8,2.0)
print(Gamma)
print(linalg.inv(Gamma))
savetxt('1e7.txt',Gamma)
savetxt('1e7_Sigma.txt',linalg.inv(Gamma))
print('1e6')
Gamma=Fisher_matrix(1e6,1e6,-0.8,1.0,0.5,3.0)
print(Gamma)
print(linalg.inv(Gamma))
savetxt('1e6.txt',Gamma)
savetxt('1e6_Sigma.txt',linalg.inv(Gamma))
print('1e5')
Gamma=Fisher_matrix(1e5,1e5,0.9,2.0,-0.8,5.0)
print(Gamma)
print(linalg.inv(Gamma))
savetxt('1e5.txt',Gamma)
savetxt('1e5_Sigma.txt',linalg.inv(Gamma))
print('1e4')
Gamma=Fisher_matrix(1e4,1e4,-0.1,3.0,-0.9,6.0)
print(Gamma)
print(linalg.inv(Gamma))
savetxt('1e4.txt',Gamma)
savetxt('1e4_Sigma.txt',linalg.inv(Gamma))

from numpy import *
from numpy.matlib import empty as empty_matrix
from scipy.integrate import *
from scipy.misc import derivative

c = 299792458
R = 149597870700
T = sqrt(R**3/(1.3271244e+20/(4*pi**2)))
varphi_0 = 0
ln_A = 0

def func_A_alpha(A_p, A_t, F_p_alpha, F_t_alpha):
    return sqrt((A_p*F_p_alpha)**2+(A_t*F_t_alpha)**2)

def func_chi_alpha(t, f_0, varphi_p_alpha, varphi_D):
    return 2*pi*f_0*t+varphi_0+varphi_p_alpha+varphi_D

def func_varphi_p_alpha(A_p, A_t, F_p_alpha, F_t_alpha):
    return arctan(-(A_t*F_t_alpha)/(A_p*F_p_alpha))

def func_varphi_D(f_0, bar_theta_S, bar_phi, bar_phi_S):
    return 2*pi*f_0/c*R*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)

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
                        *(cos(bar_theta_S)*sin(bar_theta_L)*cos(bar_phi_S)
                         -cos(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_L)))

def func_A_p(Ln):
    return 2*exp(ln_A)*(1+Ln**2)

def func_A_t(Ln):
    return -2*exp(ln_A)*Ln

def func_F_p_I(theta_S, phi_S, psi_S):
    return (cos(2*phi_S)*cos(2*psi_S)*(1+cos(theta_S)**2)/2
           -sin(2*phi_S)*sin(2*psi_S)*cos(theta_S))

def func_F_t_I(theta_S, phi_S, psi_S):
    return (cos(2*phi_S)*sin(2*psi_S)*(1+cos(theta_S)**2)/2
           -sin(2*phi_S)*cos(2*psi_S)*cos(theta_S))

def func_F_p_II(theta_S, phi_S, psi_S):
    return (sin(2*phi_S)*cos(2*psi_S)*(1+cos(theta_S)**2)/2
           +cos(2*phi_S)*sin(2*psi_S)*cos(theta_S))

def func_F_t_II(theta_S, phi_S, psi_S):
    return (sin(2*phi_S)*sin(2*psi_S)*(1+cos(theta_S)**2)/2
           +cos(2*phi_S)*cos(2*psi_S)*cos(theta_S))

def func_theta_S(bar_theta_S, bar_phi, bar_phi_S):
    return arccos(
           ((sqrt(1)/2)*cos(bar_theta_S)
           -(sqrt(3)/2)*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)))

def func_phi_S(t, bar_theta_S, bar_phi, bar_phi_S):
    alpha_0 = 0
    return alpha_0+2*pi*t/T+arctan(
        (sin(bar_theta_S)*cos(bar_phi-bar_phi_S)-sqrt(3)*cos(bar_theta_S))
        /(2*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)))

def func_psi_S(Lz, Ln, zn, nLz):
    return arctan((Lz-Ln*zn)/nLz)

def func_bar_phi(t):
    bar_phi_0 = 0
    return bar_phi_0+2*pi*t/T

def S_n(f):
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

def preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    bar_theta_S = arccos(bar_mu_S)
    bar_theta_L = arccos(bar_mu_L)
    bar_phi = func_bar_phi(t)
    theta_S = func_theta_S(bar_theta_S, bar_phi, bar_phi_S)
    Lz = func_Lz(bar_theta_L, bar_phi, bar_phi_L)
    Ln = func_Ln(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S)
    zn = func_zn(theta_S)
    nLz = func_nLz(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S, bar_phi)
    phi_S = func_phi_S(t, bar_theta_S, bar_phi, bar_phi_S)
    psi_S = func_psi_S(Lz, Ln, zn, nLz)
    A_p = func_A_p(Ln)
    A_t = func_A_t(Ln)
    return (A_p, A_t, theta_S, phi_S, psi_S)

def func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    A_p, A_t, theta_S, phi_S, psi_S = args
    F_p_I = func_F_p_I(theta_S, phi_S, psi_S)
    F_t_I = func_F_t_I(theta_S, phi_S, psi_S)
    return func_A_alpha(A_p, A_t, F_p_I, F_t_I)

def func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    A_p, A_t, theta_S, phi_S, psi_S = args
    F_p_II = func_F_p_II(theta_S, phi_S, psi_S)
    F_t_II = func_F_t_II(theta_S, phi_S, psi_S)
    return func_A_alpha(A_p, A_t, F_p_II, F_t_II)

def func_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    A_p, A_t, theta_S, phi_S, psi_S = args
    bar_theta_S, bar_phi = arccos(bar_mu_S), func_bar_phi(t)
    F_p_I = func_F_p_I(theta_S, phi_S, psi_S)
    F_t_I = func_F_t_I(theta_S, phi_S, psi_S)
    varphi_p_I = func_varphi_p_alpha(A_p, A_t, F_p_I, F_t_I)
    varphi_D = func_varphi_D(f_0, bar_theta_S, bar_phi, bar_phi_S)
    return func_chi_alpha(t, f_0, varphi_p_I, varphi_D)

def func_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    A_p, A_t, theta_S, phi_S, psi_S = args
    bar_theta_S, bar_phi = arccos(bar_mu_S), func_bar_phi(t)
    F_p_II = func_F_p_II(theta_S, phi_S, psi_S)
    F_t_II = func_F_t_II(theta_S, phi_S, psi_S)
    varphi_p_II = func_varphi_p_alpha(A_p, A_t, F_p_II, F_t_II)
    varphi_D = func_varphi_D(f_0, bar_theta_S, bar_phi, bar_phi_S)
    return func_chi_alpha(t, f_0, varphi_p_II, varphi_D)

def func_h_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    A_p, A_t, theta_S, phi_S, psi_S = args
    bar_theta_S, bar_phi = arccos(bar_mu_S), func_bar_phi(t)
    F_p_I = func_F_p_I(theta_S, phi_S, psi_S)
    F_t_I = func_F_t_I(theta_S, phi_S, psi_S)
    A_I = func_A_alpha(A_p, A_t, F_p_I, F_t_I)
    varphi_p_I = func_varphi_p_alpha(A_p, A_t, F_p_I, F_t_I)
    varphi_D = func_varphi_D(f_0, bar_theta_S, bar_phi, bar_phi_S)
    chi_I = func_chi_alpha(t, f_0, varphi_p_I, varphi_D)
    return (sqrt(3)/2)*A_I*cos(chi_I)

def func_h_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    args = preparation(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    A_p, A_t, theta_S, phi_S, psi_S = args
    bar_theta_S, bar_phi = arccos(bar_mu_S), func_bar_phi(t)
    F_p_II = func_F_p_II(theta_S, phi_S, psi_S)
    F_t_II = func_F_t_II(theta_S, phi_S, psi_S)
    A_II = func_A_alpha(A_p, A_t, F_p_II, F_t_II)
    varphi_p_II = func_varphi_p_alpha(A_p, A_t, F_p_II, F_t_II)
    varphi_D = func_varphi_D(f_0, bar_theta_S, bar_phi, bar_phi_S)
    chi_II = func_chi_alpha(t, f_0, varphi_p_II, varphi_D)
    return (sqrt(3)/2)*A_II*cos(chi_II)

def partial_0_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)

def partial_1_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 0

def partial_2_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 0

def partial_3_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_I_(bar_mu_S_):
        return func_A_I(t, f_0, bar_mu_S_, bar_phi_S, bar_mu_L, bar_phi_L)
    delta = func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_I_, bar_mu_S, delta)

def partial_4_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_I_(bar_phi_S_):
        return func_A_I(t, f_0, bar_mu_S, bar_phi_S_, bar_mu_L, bar_phi_L)
    delta = func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_I_, bar_phi_S, delta)

def partial_5_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_I_(bar_mu_L_):
        return func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L_, bar_phi_L)
    delta = func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_I_, bar_mu_L, delta)

def partial_6_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_I_(bar_phi_L_):
        return func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_)
    delta = func_A_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_I_, bar_phi_L, delta)

def partial_0_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)

def partial_1_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 0

def partial_2_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 0

def partial_3_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_II_(bar_mu_S_):
        return func_A_II(t, f_0, bar_mu_S_, bar_phi_S, bar_mu_L, bar_phi_L)
    delta = func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_II_, bar_mu_S, delta)

def partial_4_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_II_(bar_phi_S_):
        return func_A_II(t, f_0, bar_mu_S, bar_phi_S_, bar_mu_L, bar_phi_L)
    delta = func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_II_, bar_phi_S, delta)

def partial_5_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_II_(bar_mu_L_):
        return func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L_, bar_phi_L)
    delta = func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_II_, bar_mu_L, delta)

def partial_6_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_A_II_(bar_phi_L_):
        return func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_)
    delta = func_A_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)*1e-8
    return derivative(func_A_II_, bar_phi_L, delta)

def partial_0_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 0

def partial_1_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 1

def partial_2_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_I_(f_0_):
        return func_chi_I(t, f_0_, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return derivative(func_chi_I_, f_0, 1e-8)

def partial_3_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_I_(bar_mu_S_):
        return func_chi_I(t, f_0, bar_mu_S_, bar_phi_S, bar_mu_L, bar_phi_L)
    return derivative(func_chi_I_, bar_mu_S, 1e-8)

def partial_4_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_I_(bar_phi_S_):
        return func_chi_I(t, f_0, bar_mu_S, bar_phi_S_, bar_mu_L, bar_phi_L)
    return derivative(func_chi_I_, bar_phi_S, 1e-8)

def partial_5_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_I_(bar_mu_L_):
        return func_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L_, bar_phi_L)
    return derivative(func_chi_I_, bar_mu_L, 1e-8)

def partial_6_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_I_(bar_phi_L_):
        return func_chi_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_)
    return derivative(func_chi_I_, bar_phi_L, 1e-8)

def partial_0_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 0

def partial_1_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    return 1

def partial_2_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_II_(f_0_):
        return func_chi_II(t, f_0_, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
    return derivative(func_chi_II_, f_0, 1e-8)

def partial_3_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_II_(bar_mu_S_):
        return func_chi_II(t, f_0, bar_mu_S_, bar_phi_S, bar_mu_L, bar_phi_L)
    return derivative(func_chi_II_, bar_mu_S, 1e-8)

def partial_4_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_II_(bar_phi_S_):
        return func_chi_II(t, f_0, bar_mu_S, bar_phi_S_, bar_mu_L, bar_phi_L)
    return derivative(func_chi_II_, bar_phi_S, 1e-8)

def partial_5_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_II_(bar_mu_L_):
        return func_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L_, bar_phi_L)
    return derivative(func_chi_II_, bar_mu_L, 1e-8)

def partial_6_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_chi_II_(bar_phi_L_):
        return func_chi_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L_)
    return derivative(func_chi_II_, bar_phi_L, 1e-8)

def signal2noise(f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    def func_integrated(t):
        h_I = func_h_I(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
        h_II = func_h_II(t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
        return h_I**2+h_II**2
    return 2*quad(func_integrated, -1/f_0, 1/f_0)[0]

def Fisher_matrix(f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L):
    Gamma = empty_matrix((7,7))
    for i in range(7):
        for j in range(7):
            def func_integrated(t):
                args = (t, f_0, bar_mu_S, bar_phi_S, bar_mu_L, bar_phi_L)
                A_I = eval(f'func_A_I(*args)')
                A_I_i = eval(f'partial_{i}_A_I(*args)')
                A_I_j = eval(f'partial_{j}_A_I(*args)')
                chi_I_i = eval(f'partial_{i}_chi_I(*args)')
                chi_I_j = eval(f'partial_{j}_chi_I(*args)')
                A_II = eval(f'func_A_II(*args)')
                A_II_i = eval(f'partial_{i}_A_II(*args)')
                A_II_j = eval(f'partial_{j}_A_II(*args)')
                chi_II_i = eval(f'partial_{i}_chi_II(*args)')
                chi_II_j = eval(f'partial_{j}_chi_II(*args)')
                return (A_I_i*A_I_j+A_I**2*(chi_I_i*chi_I_j)
                       +A_II_i*A_II_j+A_II**2*(chi_II_i*chi_II_j))
            Gamma[i,j] = (3/4)*quad(func_integrated, -1/f_0, 1/f_0)[0]
    return Gamma

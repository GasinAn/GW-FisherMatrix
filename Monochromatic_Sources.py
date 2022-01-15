from numpy import *

c = 299792458
R = 149597870700
T = sqrt(R**3/(1.3271244e+20/(4*pi**2)))
varphi_0 = 0
ln_A = 0

def A_alpha(A_p, A_t, F_p_alpha, F_t_alpha):
    return sqrt((A_p*F_p_alpha)**2+(A_t*F_t_alpha)**2)

def chi_alpha(t, f_0, varphi_p_alpha, varphi_D):
    return 2*pi*f_0*t+varphi_0+varphi_p_alpha+varphi_D

def varphi_p_alpha(A_p, A_t, F_p_alpha, F_t_alpha):
    return arctan(-(A_t*F_t_alpha)/(A_p*F_p_alpha))

def varphi_D(f_0, bar_theta_S, bar_phi, bar_phi_S):
    return 2*pi*f_0/c*R*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)

def Lz(bar_theta_L, bar_phi, bar_phi_L):
    return ((sqrt(1)/2)*cos(bar_theta_L)
           -(sqrt(3)/2)*sin(bar_theta_L)*cos(bar_phi-bar_phi_L))

def Ln(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S):
    return (cos(bar_theta_L)*cos(bar_theta_S)
           +sin(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_L-bar_phi_S))

def zn(theta_S):
    return cos(theta_S)

def nLz(bar_theta_L, bar_theta_S, bar_phi_L, bar_phi_S, bar_phi):
    return ((1/2)*sin(bar_theta_L)*sin(bar_theta_S)*sin(bar_phi_L-bar_phi_S)
            -(sqrt(3)/2)*cos(bar_phi)
                        *(cos(bar_theta_L)*sin(bar_theta_S)*sin(bar_phi_S)
                         -cos(bar_theta_S)*sin(bar_theta_L)*sin(bar_phi_L))
            -(sqrt(3)/2)*sin(bar_phi)
                        *(cos(bar_theta_S)*sin(bar_theta_L)*cos(bar_phi_S)
                         -cos(bar_theta_L)*sin(bar_theta_S)*cos(bar_phi_L)))

def A_p(Ln):
    return 2*exp(ln_A)*(1+Ln**2)

def A_t(Ln):
    return -2*exp(ln_A)*Ln

def F_p_I(theta_S, phi_S, psi_S):
    return (cos(2*phi_S)*cos(2*psi_S)*(1+cos(theta_S)**2)/2
           -sin(2*phi_S)*sin(2*psi_S)*cos(theta_S))

def F_t_I(theta_S, phi_S, psi_S):
    return (cos(2*phi_S)*sin(2*psi_S)*(1+cos(theta_S)**2)/2
           -sin(2*phi_S)*cos(2*psi_S)*cos(theta_S))

def F_p_II(theta_S, phi_S, psi_S):
    return (sin(2*phi_S)*cos(2*psi_S)*(1+cos(theta_S)**2)/2
           +cos(2*phi_S)*sin(2*psi_S)*cos(theta_S))

def F_t_II(theta_S, phi_S, psi_S):
    return (sin(2*phi_S)*sin(2*psi_S)*(1+cos(theta_S)**2)/2
           +cos(2*phi_S)*cos(2*psi_S)*cos(theta_S))

def theta_S(bar_theta_S, bar_phi, bar_phi_S):
    return arccos(
           ((sqrt(1)/2)*cos(bar_theta_S)
           -(sqrt(3)/2)*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)))

def phi_S(t, bar_theta_S, bar_phi, bar_phi_S):
    alpha_0 = 0
    return alpha_0+2*pi*t/T+arctan(
        (sin(bar_theta_S)*cos(bar_phi-bar_phi_S)-sqrt(3)*cos(bar_theta_S))
        /(2*sin(bar_theta_S)*cos(bar_phi-bar_phi_S)))

def psi_S(Lz, Ln, zn, nLz):
    return arctan((Lz-Ln*zn)/nLz)

def phi(t):
    phi_0 = 0
    return phi_0+2*pi*t/T

def S_n(f):
    alpha = 10**(-22.79)*(f/1e-3)**(-7/3)
    beta  = 10**(-24.54)*(f/1e-3)
    gamma = 10**(-23.04)
    S_n_in = 5.049e+5*(alpha*2+beta*2+gamma*2)
    if f<=10**(-3.15):
        S_n_co = 10**(-42.685)*f**(-1.9)
    elif f<=10**(-2.75):
        S_n_co = 10**(-60.325)*f**(-7.5)
    else:
        S_n_co = 10**(-46.850)*f**(-2.6)
    return S_n_in+S_n_co

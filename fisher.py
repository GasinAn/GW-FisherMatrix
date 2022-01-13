from numpy.matlib import empty as empty_matrix
from scipy.integrate import quad


def inner_product(g, h, Sn, f_min, f_max):
    '''
    calculate the inner product.
    $$ 4 \int_{f_min}^{f_max} [Re(g)Re(h)+Im(g)Im(h)]/Sn df $$
    '''

    def function_integrated(f):
        return (g(f).real*h(f).real+g(f).imag*h(f).imag)/Sn(f)

    return 4*quad(function_integrated, f_min, f_max)[0]

def Fisher_matrix(derivs, Sn, f_min, f_max):
    '''
    calculate Fisher matrix.
    '''

    n = len(derivs)
    Gamma = empty_matrix((n,n))

    for i in range(0,n-1):
        for j in range(i+1,n):
            Gamma[i,j] = inner_product(derivs[i], derivs[j], Sn, f_min, f_max)
            Gamma[j,i] = Gamma[i,j]
    for k in range(n):
        Gamma[k,k] = inner_product(derivs[k], derivs[k], Sn, f_min, f_max)

    return Gamma


if (__name__ == '__main__'):
    from math import pi, sqrt
    from astropy.constants import G, M_sun, c

    M_sun = (M_sun*G/c**3).value
    f_0 = 70

    def Poisson_Will_1995(m1, m2, table, prior, epsilon=1):
        '''
        calculate table2 & table3 in Poisson & Will, 1995.
        '''

        m1, m2 = m1*M_sun, m2*M_sun

        M = m1+m2
        mu = m1*m2/(m1+m2)
        eta = mu/M
        calM = eta**(3/5)*M

        f_min, f_max = 10, (6**(3/2)*pi*M)**(-1)

        rou, beta, sigma = 10, 0, 0

        def function_integrated(f):
            return f**(-7/3)/(f**(-4)+2+2*f**2)
        I_7 = quad(function_integrated, f_min/f_0, f_max/f_0)[0]

        calA = sqrt(rou**2/(20*f_0**(-4/3)*I_7))

        def tilde_h(f):
            return calA*f**(-7/6)

        def deriv0(f):
            return tilde_h(f)

        def deriv1(f):
            return (2*pi*1j)*(f/f_0)*tilde_h(f)

        def deriv2(f):
            return (-1j)*tilde_h(f)

        def deriv3(f):
            v = (pi*M*f)**(1/3)
            A_4 = (
                (4/3)
                *((743/336)+(11/4)*eta)
            )
            B_4 = (
                (8/5)
                *((4*pi)-beta)
            )
            C_4 = (
                (2*epsilon)
                *((3058673/1016064)+(5429/1008)*eta+(617/144)*eta**2-sigma)
            )
            return (
                -(5j/128)*(pi*calM*f)**(-5/3)
                *(1+A_4*v**2-B_4*v**3+C_4*v**4)*tilde_h(f)
            )

        def deriv4(f):
            v = (pi*M*f)**(1/3)
            A_5 = (
                (1)
                *((743/168)-(33/4)*eta)
            )
            B_5 = (
                (27/5)
                *((4*pi)-beta)
            )
            C_5 = (
                (18*epsilon)
                *((3058673/1016064)-(5429/4032)*eta-(617/96)*eta**2-sigma)
            )
            return (
                -(1j/96)*(pi*calM*f)**(-5/3)
                *(A_5*v**2-B_5*v**3+C_5*v**4)*tilde_h(f)
            )

        def deriv5(f):
            return (3j/32)*eta**(-3/5)*(pi*calM*f)**(-2/3)*tilde_h(f)

        def deriv6(f):
            return -(15j/64)*eta**(-4/5)*(pi*calM*f)**(-1/3)*tilde_h(f)

        def Sn(f):
            return (1/5)*((f_0/f)**4+2+2*(f/f_0)**2)

        derivs = (deriv0, deriv1, deriv2, deriv3, deriv4, deriv5, deriv6)
        Gamma = Fisher_matrix(derivs, Sn, f_min, f_max)

        if prior:
            Gamma[5,5] = Gamma[5,5]+(1/8.5)**2
            Gamma[6,6] = Gamma[6,6]+(1/5.0)**2

        if (table == 2):
            Sigma = Gamma[:6,:6].I
            estimation = {}
            estimation['prior'] = prior
            estimation['epsilon'] = epsilon
            estimation['Dt_c(ms)'] = sqrt(Sigma[1,1])/f_0*1000
            estimation['Dphi_c'] = sqrt(Sigma[2,2])
            estimation['DcalM/calM(%)'] = sqrt(Sigma[3,3])*100
            estimation['Deta/eta'] = sqrt(Sigma[4,4])
            estimation['Dbeta'] = sqrt(Sigma[5,5])
            estimation['c_calMeta'] = Sigma[3,4]/sqrt(Sigma[3,3]*Sigma[4,4])
            estimation['c_calMbeta'] = Sigma[3,5]/sqrt(Sigma[3,3]*Sigma[5,5])
            estimation['c_etabeta'] = Sigma[4,5]/sqrt(Sigma[4,4]*Sigma[5,5])
            return estimation

        if (table == 3):
            Sigma = Gamma[:7,:7].I
            estimation = {}
            estimation['prior'] = prior
            estimation['Dt_c(ms)'] = sqrt(Sigma[1,1])/f_0*1000
            estimation['Dphi_c'] = sqrt(Sigma[2,2])
            estimation['DcalM/calM(%)'] = sqrt(Sigma[3,3])*100
            estimation['Deta/eta'] = sqrt(Sigma[4,4])
            estimation['Dbeta'] = sqrt(Sigma[5,5])
            estimation['Dsigma'] = sqrt(Sigma[6,6])
            estimation['c_calMeta'] = Sigma[3,4]/sqrt(Sigma[3,3]*Sigma[4,4])
            estimation['c_calMbeta'] = Sigma[3,5]/sqrt(Sigma[3,3]*Sigma[5,5])
            estimation['c_calMsigma'] = Sigma[3,6]/sqrt(Sigma[3,3]*Sigma[6,6])
            estimation['c_etabeta'] = Sigma[4,5]/sqrt(Sigma[4,4]*Sigma[5,5])
            estimation['c_etasigma'] = Sigma[4,6]/sqrt(Sigma[4,4]*Sigma[6,6])
            estimation['c_betasigma'] = Sigma[5,6]/sqrt(Sigma[5,5]*Sigma[6,6])
            return estimation

    print('Table2')
    print('SystemA')
    print(Poisson_Will_1995(m1=1.4, m2=1.4, table=2, prior=True, epsilon=1))
    print(Poisson_Will_1995(m1=1.4, m2=1.4, table=2, prior=False, epsilon=1))
    print(Poisson_Will_1995(m1=1.4, m2=1.4, table=2, prior=False, epsilon=0))
    print(Poisson_Will_1995(m1=1.4, m2=1.4, table=2, prior=False, epsilon=-1))
    print('SystemB')
    print(Poisson_Will_1995(m1=1.4, m2=10, table=2, prior=True, epsilon=1))
    print(Poisson_Will_1995(m1=1.4, m2=10, table=2, prior=False, epsilon=1))
    print(Poisson_Will_1995(m1=1.4, m2=10, table=2, prior=False, epsilon=0))
    print(Poisson_Will_1995(m1=1.4, m2=10, table=2, prior=False, epsilon=-1))
    print('SystemC')
    print(Poisson_Will_1995(m1=10, m2=10, table=2, prior=True, epsilon=1))
    print(Poisson_Will_1995(m1=10, m2=10, table=2, prior=False, epsilon=1))
    print(Poisson_Will_1995(m1=10, m2=10, table=2, prior=False, epsilon=0))
    print(Poisson_Will_1995(m1=10, m2=10, table=2, prior=False, epsilon=-1))

    print('Table3')
    print('SystemA')
    print(Poisson_Will_1995(m1=1.4, m2=1.4, table=3, prior=True))
    print(Poisson_Will_1995(m1=1.4, m2=1.4, table=3, prior=False))
    print('SystemB')
    print(Poisson_Will_1995(m1=1.4, m2=10, table=3, prior=True))
    print(Poisson_Will_1995(m1=1.4, m2=10, table=3, prior=False))
    print('SystemC')
    print(Poisson_Will_1995(m1=10, m2=10, table=3, prior=True))
    print(Poisson_Will_1995(m1=10, m2=10, table=3, prior=False))

import numpy as np


for sys in ('A', 'B', 'C'):

    Gamma = np.loadtxt(f'sys{sys}_epsilon=1_Gamma.txt')[:6,:6]
    Gamma[5,5] = Gamma[5,5]+(1/8.5)**2 
    Sigma = np.linalg.inv(Gamma)
    s_2 = str(np.sqrt(Sigma[1,1])/70*1000)
    s_3 = str(np.sqrt(Sigma[2,2]))
    s_4 = str(np.sqrt(Sigma[3,3])*100)
    s_5 = str(np.sqrt(Sigma[4,4]))
    s_6 = str(np.sqrt(Sigma[5,5]))
    c_45 = str(Sigma[3,4]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[4,4])))
    c_46 = str(Sigma[3,5]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[5,5])))
    c_56 = str(Sigma[4,5]/(np.sqrt(Sigma[4,4])*np.sqrt(Sigma[5,5])))
    print(f'{s_2} {s_3} {s_4} {s_5} {s_6} {c_45} {c_46} {c_56}\n')

    Gamma = np.loadtxt(f'sys{sys}_epsilon=1_Gamma.txt')[:6,:6]
    Sigma = np.linalg.inv(Gamma)
    s_2 = str(np.sqrt(Sigma[1,1])/70*1000)
    s_3 = str(np.sqrt(Sigma[2,2]))
    s_4 = str(np.sqrt(Sigma[3,3])*100)
    s_5 = str(np.sqrt(Sigma[4,4]))
    s_6 = str(np.sqrt(Sigma[5,5]))
    c_45 = str(Sigma[3,4]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[4,4])))
    c_46 = str(Sigma[3,5]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[5,5])))
    c_56 = str(Sigma[4,5]/(np.sqrt(Sigma[4,4])*np.sqrt(Sigma[5,5])))
    print(f'{s_2} {s_3} {s_4} {s_5} {s_6} {c_45} {c_46} {c_56}\n')

    Gamma = np.loadtxt(f'sys{sys}_epsilon=0_Gamma.txt')[:6,:6]
    Sigma = np.linalg.inv(Gamma)
    s_2 = str(np.sqrt(Sigma[1,1])/70*1000)
    s_3 = str(np.sqrt(Sigma[2,2]))
    s_4 = str(np.sqrt(Sigma[3,3])*100)
    s_5 = str(np.sqrt(Sigma[4,4]))
    s_6 = str(np.sqrt(Sigma[5,5]))
    c_45 = str(Sigma[3,4]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[4,4])))
    c_46 = str(Sigma[3,5]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[5,5])))
    c_56 = str(Sigma[4,5]/(np.sqrt(Sigma[4,4])*np.sqrt(Sigma[5,5])))
    print(f'{s_2} {s_3} {s_4} {s_5} {s_6} {c_45} {c_46} {c_56}\n')

    Gamma = np.loadtxt(f'sys{sys}_epsilon=-1_Gamma.txt')[:6,:6]
    Sigma = np.linalg.inv(Gamma)
    s_2 = str(np.sqrt(Sigma[1,1])/70*1000)
    s_3 = str(np.sqrt(Sigma[2,2]))
    s_4 = str(np.sqrt(Sigma[3,3])*100)
    s_5 = str(np.sqrt(Sigma[4,4]))
    s_6 = str(np.sqrt(Sigma[5,5]))
    c_45 = str(Sigma[3,4]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[4,4])))
    c_46 = str(Sigma[3,5]/(np.sqrt(Sigma[3,3])*np.sqrt(Sigma[5,5])))
    c_56 = str(Sigma[4,5]/(np.sqrt(Sigma[4,4])*np.sqrt(Sigma[5,5])))
    print(f'{s_2} {s_3} {s_4} {s_5} {s_6} {c_45} {c_46} {c_56}\n')

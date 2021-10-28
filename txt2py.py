for sys in ('A', 'B', 'C'):
    for epsilon in ('1', '0', '-1'):
        txt_file_name = f'sys{sys}_epsilon={epsilon}.txt'
        with open(txt_file_name, 'r', encoding='utf-16') as txt_file:
            txt = txt_file.readlines()
            with open(f'sys{sys}_epsilon={epsilon}.py', 'w') as py_file:
                py_file.write('from math import pi\n'
                              '\n'
                              'import numpy as np\n'
                              'from astropy.constants import G, M_sun, c\n'
                              'from scipy.integrate import quad\n'
                              '\n'
                )
                py_file.write('\n')
                n = 1
                for i in range(1,8):
                    for j in range(i,8):
                        py_file.write(f'def prod_{i}{j}(f):\n')
                        py_file.write(f'    return {txt[4*n-3]}\n')
                        n = n+1
                py_file.write('\n')

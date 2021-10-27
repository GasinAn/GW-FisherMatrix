for sys in ('A', 'B', 'C'):
    for epsilon in ('1', '0', '-1'):
        with open(f'sys{sys}_epsilon={epsilon}.txt', 'r') as txt_file:
            txt = txt_file.readlines()
        with open(f'sys{sys}_epsilon={epsilon}.txt', 'w') as txt_file:
            for txt_line in txt:
                if txt_line != '\n':
                    txt_file.write(txt_line)
        with open(f'sys{sys}_epsilon={epsilon}.txt', 'r') as txt_file:
            txt = txt_file.readlines()
            with open(f'sys{sys}_epsilon={epsilon}.py', 'w') as py_file:
                py_file.write('from math import pi\n\n\n')
                n = 1
                for i in range(1,8):
                    for j in range(i,8):
                        py_file.write(f'def prod_{i}{j}(f):\n')
                        py_file.write(f'    return {txt[2*n-1]}\n')
                        n = n+1

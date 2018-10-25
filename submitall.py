#! /usr/bin/env python

import sys, os
import subprocess
from time import sleep

if sys.argv[1]:
    my_path = sys.argv[1]
else:
    my_path = os.getcwd()

for root, dirs, files in os.walk(my_path):
#    with open("/projects/bioenergy/mjankous/Mo2C_6_clean/results/adsorbates/H/clean/INCAR",'r') as h:
#    root_incar = h.read()
#    with open("/projects/bioenergy/mjankous/Mo2C_6_clean/results/adsorbates/H/clean/POTCAR",'r') as j:
#    root_potcar = j.read()
#    with open("/projects/bioenergy/mjankous/Mo2C_6_clean/results/adsorbates/H/clean/KPOINTS",'r') as k:
#    root_kpts = k.read()
    for name in dirs:
        output_path = os.path.join(root, name)
        if os.path.exists(output_path + '/POSCAR~'):
            with open(output_path + '/POSCAR~','r') as f:
                orig_poscar = f.read()
            with open(output_path + '/POSCAR','w') as g:
                g.write(orig_poscar)

#    if os.path.exists(paff + '/INCAR'):
#        os.rename(paff + '/INCAR', paff + '/old_INCAR')
#        with open(paff + '/old_INCAR','r') as f:
#            my_incar = f.readlines()
#            new_incar_lst = []
#            for a_line in my_incar:
#                if a_line.strip():
#                    lin_data = a_line.split()
#                        if lin_data[0] == 'LREAL':
#                            print('LREAL FOUND')
#                            lin_data = ['LREAL', '=', 'Auto']
#                    new_line = ' '.join(lin_data)
#                    new_incar_lst.append(new_line)
#        new_incar = '\n'.join(new_incar_lst)
#        with open(paff + '/INCAR','w') as g:
#            g.write(new_incar)
            
        if os.path.exists(output_path + '/POSCAR'):
#        with open(paff + '/POTCAR','w') as pot_to_write:
#            pot_to_write.write(root_potcar)
#        with open(paff + '/INCAR','w') as in_to_write:
#            in_to_write.write(root_incar)
#        with open(paff + '/KPOINTS','w') as kpts_to_write:
#            kpts_to_write.write(root_kpts)
            jobname_lst = output_path.split('/')
            jobname = '_'.join(jobname_lst[6:])
            jobname += '.log'
            print(jobname)
            os.chdir(output_path)
#        print(os.getcwd())
            subprocess.Popen(['qbio_e', '1', '24', '24', 'batch', jobname])
#qbio_e is my personal submission script, 1 refers to number of nodes, first 24 is number of cores, second 24 is number of hours, batch is the queue
            sleep(2)

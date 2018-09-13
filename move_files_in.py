#! /usr/bin/env python

import sys, os

if sys.argv[1]:
    my_path = sys.argv[1]
else:
    my_path = os.getcwd()

with open(my_path + '/POTCAR','r') as potcar_origin:
    potcar = potcar_origin.read()

with open(my_path + '/INCAR','r') as incar_origin:
    incar = incar_origin.read()

with open(my_path + '/KPOINTS','r') as kpts_origin:
    kpts = kpts_origin.read()

for root, dirs, files in os.walk(my_path):
    for name in dirs:
        output_path = os.path.join(root, name)
        if os.path.exists(output_path + '/POSCAR'):
            with open(output_path + '/POTCAR','w') as potcar_destination:
                potcar_destination.write(potcar)
            with open(output_path + '/INCAR','w') as incar_destination:
                incar_destination.write(incar)
            with open(output_path + '/KPOINTS','w') as kpts_destination:
                kpts_destination.write(kpts)

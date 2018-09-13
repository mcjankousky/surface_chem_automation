#! /usr/bin/env python

import sys, os

def move_files_in(top_dir, calc_dir):
    
    with open(top_dir + '/POTCAR','r') as potcar_origin:
        potcar = potcar_origin.read()
    
    with open(top_dir + '/INCAR','r') as incar_origin:
        incar = incar_origin.read()
    
    with open(top_dir + '/KPOINTS','r') as kpts_origin:
        kpts = kpts_origin.read()
        
    walk_dir=top_dir+calc_dir
    
    for root, dirs, files in os.walk(walk_dir):
        for name in dirs:
            output_path = os.path.join(root, name)
            if os.path.exists(output_path + '/POSCAR'):
                with open(output_path + '/POTCAR','w') as potcar_destination:
                    potcar_destination.write(potcar)
                with open(output_path + '/INCAR','w') as incar_destination:
                    incar_destination.write(incar)
                with open(output_path + '/KPOINTS','w') as kpts_destination:
                    kpts_destination.write(kpts)
    return
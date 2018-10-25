# /usr/bin/env python

import sys, os
import pandas as pd
import numpy as np

def E_conv_collect(calc_dir,calc_name,surf_stability=True):
    calc_lst = []
    E_lst = []
    if surf_stability:
        species_lst = []
        area_lst = []
    for root, dirs, files in os.walk(calc_dir):
        for name in dirs:
            calc_path = os.path.join(root, name)
            if os.path.exists(calc_path + '/OUTCAR'):
                calc_lst.append(calc_path)
                conv = False
                with open(calc_path + '/OUTCAR') as outcar:
                    if 'Elaps' in outcar.read():
                        conv = True
                #conv_lst.append(conv)	
                if conv == True:
                    if os.path.exists(calc_path + '/OSZICAR'):
                        with open(calc_path + '/OSZICAR') as oszicar:
                            oszicar = oszicar.readlines()
                            if len(oszicar) > 0:
                                last_line = oszicar[-1]
                                last_line = last_line.split()
                                last_line = [x for x in last_line if x != '']
                                Energy = last_line[2]
                            else:
                                Energy = 0
                    else:
                        Energy = 0
                else:
                    Energy = 0
                E_lst.append(Energy)
                if surf_stability:
                    if os.path.exists(calc_path + '/POSCAR'):
                        with open(calc_path + '/POSCAR') as poscar:
                            poscar_lines = poscar.readlines()
                            species_str = poscar_lines[5]
                            species_num = poscar_lines[6]
                            species_types = species_str.split()
                            species_types = [x for x in species_types if x != '']
                            species_num_lst = species_num.split()
                            species_num_lst = [x for x in species_num_lst if x != '']
                            species_dict = {}
                            for idx in range(0, len(species_types)):
                                species_dict[species_types[idx]] = int(species_num_lst[idx])
                            a_vector = poscar_lines[2]
                            b_vector = poscar_lines[3]
                            a_vector = a_vector.split()
                            a_vector = [float(x) for x in a_vector if x != '']
                            b_vector = b_vector.split()
                            b_vector = [float(x) for x in b_vector if x != '']
                            area = np.linalg.norm(np.cross(a_vector, b_vector))
                        species_lst.append(species_dict)
                        area_lst.append(area)
                        
                    else:
                        species_lst.append('no POSCAR found')
                            
                    
    				
    E_dict = {}
    E_dict['calculation'] = calc_lst
    E_dict['E'] = E_lst
    
    if surf_stability:
        E_dict['species'] = species_lst
        E_dict['area'] = area_lst
    
    E_df = pd.DataFrame.from_dict(E_dict)
    
    if surf_stability:
        E_df = E_df[['calculation','E','species','area']]
    else:
        E_df = E_df[['calculation','E']]
    
    E_df.sort_values(by='E', ascending=True, inplace=True) # DVF comment: Sorting so we can read-in lowest energy structures in next step by taking top rows
    
    E_df.to_csv(calc_name + '.csv')

    return E_df
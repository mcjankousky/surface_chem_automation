# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 14:16:49 2018

@author: mjankous
"""
import os
import pymatgen as mg
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.surface import get_symmetrically_distinct_miller_indices
from pymatgen.analysis.adsorption import AdsorbateSiteFinder

"""
The script below provides a framework to create all of the surfaces with a given
maximum miller index. It is very similer to the existing generate_all_slabs 
functionality within pymatgen with a few changes. First, it uses my 'equivalent_surfaces'
symmetry and "make_single_species_termination" method to try to enumerate those
surfaces previously missed by pymatgen
Second, it makes the slabs with varying thicknesses so that convergence of 
surface energies with respect to layers can be tested. 
Third, it makes surfaces as 2x2x1 supercells by default, as carrie suggested the 
extra room might give us better convergence as the surfaces relax
Fourth, it freezes all of the atoms more than 5 angstroms from the surface of the 
slab to help the calculations run more quickly.
Finally, it also includes a framework to save all of the generated surfaces.
"""

mat_path='\\Mo2C\\orthorhombic'

struc = Poscar.from_file(os.getcwd() + '%s\\bulk\\CONTCAR' %mat_path).structure

#miller = [0,0,1]

def freeze_center(slab):
    sf = AdsorbateSiteFinder(slab)
    surf_sites = sf.find_surface_sites_by_height(slab, height=5, bottom=True)
    sd_lst = []
    for site in slab:
        if site in surf_sites:
            sd_lst.append([True, True, True])
        else:
            sd_lst.append([False, False, False])
    slab.add_site_property('selective_dynamics',sd_lst)
            
    return slab
            
slab_memory = []
            
miller_lst = get_symmetrically_distinct_miller_indices(struc,1)
for miller in miller_lst:
    for n_angstroms in range(22, 28):
        slabgen = SlabGenerator(struc,miller,n_angstroms,15, in_unit_planes=False, lll_reduce=False)

        slabs = slabgen.get_slabs(symmetrize='equivalent_surface')

        slab_lst = []

        for slab in slabs:
                make_term = slab.make_single_species_termination()
                slab_lst.append(slab)
                if type(make_term) == list:
                    slab_lst.extend(make_term)
        
        i = 0

        for slab in slab_lst:
            new_slab = slab.copy()
            new_slab = freeze_center(new_slab)
            new_slab.get_sorted_structure()
            if new_slab in slab_memory:
                pass
            else:                
                i += 1
                if not os.path.exists(os.getcwd() + '%s\surface_termination\%s%s%s\%s\%sangstroms' %(mat_path,miller[0],miller[1],miller[2],i,n_angstroms)):
                    os.makedirs(os.getcwd() + '%s\surface_termination\%s%s%s\%s\%sangstroms' %(mat_path,miller[0],miller[1],miller[2],i,n_angstroms))
                new_slab.to('poscar',os.getcwd() + '%s\surface_termination\%s%s%s\%s\%sangstroms\\POSCAR' %(mat_path,miller[0],miller[1],miller[2],i,n_angstroms))
            slab_memory.append(new_slab)
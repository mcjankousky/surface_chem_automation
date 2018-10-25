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
Third, it freezes all of the atoms more than 5 angstroms from the surface of the 
slab to help the calculations run more quickly.
Finally, it also includes a framework to save all of the generated surfaces.
"""

#mat_path='\\Mo2C\\orthorhombic'
#
#struc = Poscar.from_file(os.getcwd() + '%s\\bulk\\CONTCAR' %mat_path).structure

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

def get_all_slabs(unit_cell, max_miller_ind, slab_thickness, surf_supercell, run_dir, in_unit_planes=False):
    slab_memory = []
    if not in_unit_planes:
        slab_range = range(slab_thickness-3, slab_thickness+3)
        layer_units = 'angstroms'
    else:
        slab_range = range(slab_thickness-1,slab_thickness+1)
        layer_units = 'layers'
    miller_lst = get_symmetrically_distinct_miller_indices(unit_cell, max_miller_ind)
    for miller in miller_lst:
        for n_angstroms in slab_range:
            slabgen = SlabGenerator(unit_cell,miller,n_angstroms,15, in_unit_planes=in_unit_planes, lll_reduce=True)
    
            slabs = slabgen.get_slabs(symmetrize='equivalent_surface')
    
            slab_lst = []
    
            for slab in slabs:
                    make_term = slab.make_single_species_termination()
                    slab_lst.append(slab)
                    if type(make_term) == list:
                        slab_lst.extend(make_term)
            
            slab_idx = 0
    
            for slab in slab_lst:
                new_slab = slab.copy()
                new_slab = freeze_center(new_slab)
                new_slab.make_supercell([surf_supercell[0],surf_supercell[1],1])
                new_slab = new_slab.get_sorted_structure()
                if new_slab in slab_memory:
                    pass
                else:                
                    slab_idx += 1
                    if not os.path.exists('%s\surface_stability\%s%s%s\%s\%s%s' %(run_dir,miller[0],miller[1],miller[2],slab_idx,n_angstroms,layer_units)):
                        os.makedirs('%s\surface_stability\%s%s%s\%s\%s%s' %(run_dir,miller[0],miller[1],miller[2],slab_idx,n_angstroms,layer_units))
                    new_slab.to('poscar','%s\surface_stability\%s%s%s\%s\%s%s\\POSCAR' %(run_dir,miller[0],miller[1],miller[2],slab_idx,n_angstroms,layer_units))
                    slab_memory.append(new_slab)
    if not os.path.exists('%s/surface_stability/bulk_reference' %run_dir):
        os.makedirs('%s/surface_stability/bulk_reference' %run_dir)
    unit_cell.to("poscar",'%s/surface_stability/bulk_reference/POSCAR' %run_dir)
    return slab_memory
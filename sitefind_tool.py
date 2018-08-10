# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 16:29:57 2018

@author: mjankous
"""

import pymatgen as mg
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.surface import SlabGenerator
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
import itertools
import os

"""
The script below provides a framework which can be used to create slabs with
varying coverage conditions from a root slab that has been generated elsewhere
"""


slab = Poscar.from_file(os.getcwd() + "\\Mo2C\\orthorhombic\\example_slab\\POSCAR").structure
#a choice of root slab that I have been using to test out the framework below



def generate_site_lst(slab):
    sf = AdsorbateSiteFinder(slab)
    #creates an AdsorbateSiteFinder object from pymatgen.analysis.adsorption 
    #to identify the possible sites for adsorbates    
    
    dict_ads_site = sf.find_adsorption_sites(distance=1.7, symm_reduce=None)
    #a dictionary with all possible adsorption sites on the surface
    
    site_lst = []
    
    for site_type in dict_ads_site.keys():
        i = 0
        if site_type != 'all':
            for site in dict_ads_site[site_type]:
                i += 1
                site_name = site_type + '_' + str(i)
                site_lst.append(site_name)
    #a list of all the sites with somewhat intuitive names, like hollow_1, ontop_2, or bridge_3
    return site_lst, dict_ads_site

def make_CO(cdown=True):
    #Creates a carbon monoxide molecule to adsorb to the surface. 
    #Always oriented normal to the surface, if cdown is True the C atom is closer to the surface
    if cdown == True:
        mol = mg.core.structure.Molecule(['C','O'],[np.array([0, 0, 0]),np.array([0, 0, 1.13])])
    else:
        mol = mg.core.structure.Molecule(['O','C'],[np.array([0, 0, 0]), np.array([0, 0, 1.13])])
    return mol    
    
def make_NH3():
    #creates an ammonia molecule to adsorb to the surface
    #only one orientation, because I've mostly only seen this configuration for NH3 on surfaces
    mol = mg.core.structure.Molecule(['N','H','H','H'], [np.array([0, 0, 0]),
                                                         np.array([-0.8248, -0.4762, 0.336]),
                                                         np.array([0.8248, -0.4762, 0.336]),
                                                         np.array([0, 0.9524, 0.336])])
    return mol

def make_H():
    #creates an ammonia molecule to adsorb to the surface
    #only one orientation, because I've mostly only seen this configuration for NH3 on surfaces
    mol = mg.core.structure.Molecule(['H'],[np.array([0,0,0])])
    return mol

def make_site_combos(slab,n_sites):
    #takes all of the possible adsorption site labels and makes all possible
    #combinations of a number n_sites of the sites with no repetition
    site_lst, dict_ads_site= generate_site_lst(slab)
    site_combos = []
    for combo in itertools.combinations(site_lst,n_sites):
        site_combos.append(combo)
    return site_combos


def save_site_combos(slab, adsorbate, path, n_sites, dist_reduce=2.1):
    site_combos = make_site_combos(slab,n_sites)
    #generates combinations of sites
    site_lst, dict_ads_site = generate_site_lst(slab)
    #generates individual sites
    idx_lst = []
    for idx in range(1, n_sites+1):
        neg_idx = idx * (-1)
        idx_lst.append(neg_idx)
    #a list of the indices for adsorbates used for calculating distance between adsorbates
    for combo in site_combos:
        dirname_lst = []
        init_slab = slab.copy()
        selected_sites = []
        for sitename in combo:
            dirname_lst.append(sitename)
            sitename_splt = sitename.split('_')
            site_type = sitename_splt[0]
            site_id = int(sitename_splt[1]) - 1
            site = dict_ads_site[site_type][site_id]
            init_slab.append('H', site, coords_are_cartesian=True)
            selected_sites.append(site)
            #stores a list of sites for a combination and appends hydrogen
            #in those sites to perform distance measurements between the sites
        dirname = '_'.join(dirname_lst)
        dist_lst = []
        for idx_combo in itertools.combinations(idx_lst,2):
            dist = init_slab.get_distance(idx_combo[0], idx_combo[1])
            dist_lst.append(dist)
            #measures the distance between the sites stored above
        if (n_sites == 1) or (dist_reduce < min(dist_lst)):
            fin_slab = slab.copy()
            for site in selected_sites:
                sf = AdsorbateSiteFinder(fin_slab)
                fin_slab = sf.add_adsorbate(adsorbate,site)
                fin_slab = fin_slab.get_sorted_structure()
                #appends the specified adsorbate to the slab in the selected 
                #sites if the distance between the sites is more than a 
                #specified number of angstroms
            if not os.path.exists(os.getcwd() + '\\Mo2C\\orthorhombic\\example_slab\\adsorbates\\%s%s'  %(path, dirname)):
                os.makedirs(os.getcwd() + '\\Mo2C\\orthorhombic\\example_slab\\adsorbates\\%s%s' %(path, dirname))
            fin_slab.to('poscar',os.getcwd() + '\\Mo2C\\orthorhombic\\example_slab\\adsorbates\\%s%s\\POSCAR' %(path, dirname))
            #stores the generated slabs with adsorbates in a local directory


def halfMLproximity(slab):
    prox_dict = {}
    radii = [(0,0.5),(0.5,0.7),(0.7,0.9),(0.9,1.1),(1.1,1.3),(1.3,1.5),(1.5,1.7),
             (1.7,1.9),(1.9,2.1),(2.1,2.3),(2.3,2.5),(2.5,2.7),(2.7,2.9),(2.9,3.1),
             (3.1,3.3), (3.3,3.5),(3.5,3.7),(3.7,3.9),(3.9,10)]
    halfMLcombos = make_site_combos(slab,2)
    site_lst, dict_ads_site = generate_site_lst(slab)
    for radius in radii:
        case_lst = []
        for combo in halfMLcombos:
            dirname = combo[0] + '_' + combo[1]
            sitename_splt1 = combo[0].split('_')
            sitename_splt2 = combo[1].split('_')
            site_type1 = sitename_splt1[0]
            site_id1 = int(sitename_splt1[1]) - 1
            site_type2 = sitename_splt2[0]
            site_id2 = int(sitename_splt2[1]) - 1
            site1 = dict_ads_site[site_type1][site_id1]
            site2 = dict_ads_site[site_type2][site_id2]
            slab = Poscar.from_file(os.getcwd() + "\\Mo2C\\orthorhombic\\example_slab\\POSCAR").structure
            slab.append('H',site1,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            slab.append('H',site2,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            if radius[0] <= slab.get_distance(len(slab) - 2, len(slab) - 1) <= radius[1]:
                case_lst.append(dirname)
        prox_dict[str(radius)] = case_lst
    return prox_dict

def threequarterMLproximity(slab):
    prox_dict = {}
    radii = [(0,0.5),(0.5,0.7),(0.7,0.9),(0.9,1.1),(1.1,1.3),(1.3,1.5),(1.5,1.7),
             (1.7,1.9),(1.9,2.1),(2.1,2.3),(2.3,2.5),(2.5,2.7),(2.7,2.9),(2.9,3.1),
             (3.1,3.3), (3.3,3.5),(3.5,3.7),(3.7,3.9),(3.9,10)]
    threequarterMLcombos = make_site_combos(slab,3)
    site_lst, dict_ads_site = generate_site_lst(slab)
    for radius in radii:
        case_lst = []
        k = 0
        for combo in threequarterMLcombos:
            dirname = combo[0] + '_' + combo[1] + '_' + combo[2]
            while k < 5:
                print(dirname)
                k += 1
            sitename_splt1 = combo[0].split('_')
            sitename_splt2 = combo[1].split('_')
            sitename_splt3 = combo[2].split('_')
            site_type1 = sitename_splt1[0]
            site_id1 = int(sitename_splt1[1]) - 1
            site_type2 = sitename_splt2[0]
            site_id2 = int(sitename_splt2[1]) - 1
            site_type3 = sitename_splt3[0]
            site_id3 = int(sitename_splt3[1]) - 1
            site1 = dict_ads_site[site_type1][site_id1]
            site2 = dict_ads_site[site_type2][site_id2]
            site3 = dict_ads_site[site_type3][site_id3]
            slab = Poscar.from_file(os.getcwd() + "\\Mo2C\\orthorhombic\\ex_slab\\POSCAR").structure
            slab.append('H',site1,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            slab.append('H',site2,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            slab.append('H',site3,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            ads_distances = [slab.get_distance(len(slab) - 2, len(slab) - 1),
                             slab.get_distance(len(slab) - 2, len(slab) - 3),
                             slab.get_distance(len(slab) - 3, len(slab) - 1)]
            if radius[0] <= min(ads_distances) <= radius[1]:
                case_lst.append(dirname)
        prox_dict[str(radius)] = case_lst
    return prox_dict


def MLproximity(slab):
    prox_dict = {}
    radii = [(0,0.5),(0.5,0.7),(0.7,0.9),(0.9,1.1),(1.1,1.3),(1.3,1.5),(1.5,1.7),
             (1.7,1.9),(1.9,2.1),(2.1,2.3),(2.3,2.5),(2.5,2.7),(2.7,2.9),(2.9,3.1),
             (3.1,3.3), (3.3,3.5),(3.5,3.7),(3.7,3.9),(3.9,10)]
    MLcombos = make_site_combos(slab,4)
    site_lst, ads_sites = generate_site_lst(slab)
    for radius in radii:
        case_lst = []
        for combo in MLcombos:
            dirname = combo[0] + '_' + combo[1] + '_' + combo[2] + '_' + combo[3]
            sitename_splt1 = combo[0].split('_')
            sitename_splt2 = combo[1].split('_')
            sitename_splt3 = combo[2].split('_')
            sitename_splt4 = combo[3].split('_')
            site_type1 = sitename_splt1[0]
            site_id1 = int(sitename_splt1[1]) - 1
            site_type2 = sitename_splt2[0]
            site_id2 = int(sitename_splt2[1]) - 1
            site_type3 = sitename_splt3[0]
            site_id3 = int(sitename_splt3[1]) - 1
            site_type4 = sitename_splt4[0]
            site_id4 = int(sitename_splt4[1]) - 1
            site1 = ads_sites[site_type1][site_id1]
            site2 = ads_sites[site_type2][site_id2]
            site3 = ads_sites[site_type3][site_id3]
            site4 = ads_sites[site_type4][site_id4]
            slab = Poscar.from_file(os.getcwd() + "\VASP_files\Mo2C\\orthorhombic\\ex_slab\\POSCAR").structure
            slab.append('H',site1,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            slab.append('H',site2,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            slab.append('H',site3,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            slab.append('H',site4,coords_are_cartesian=True,properties={'selective_dynamics':[True, True, True]})
            ads_distances = [slab.get_distance(len(slab) - 2, len(slab) - 1),
                             slab.get_distance(len(slab) - 2, len(slab) - 3),
                             slab.get_distance(len(slab) - 3, len(slab) - 1),
                             slab.get_distance(len(slab) - 4, len(slab) - 1),
                             slab.get_distance(len(slab) - 4, len(slab) - 2),
                             slab.get_distance(len(slab) - 4, len(slab) - 3)]
            if radius[0] <= min(ads_distances) <= radius[1]:
                case_lst.append(dirname)
        prox_dict[str(radius)] = case_lst
    return prox_dict





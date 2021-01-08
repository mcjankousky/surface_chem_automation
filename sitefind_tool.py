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
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.structure import Molecule
import itertools
import os
import pandas as pd
#import networkx as nx
#import json
from scipy.spatial import distance

"""
The script below provides a framework which can be used to create slabs with
varying coverage conditions from a root slab that has been generated elsewhere
"""


#slab = Poscar.from_file("C:\\Users\mjankous\Documents\VASP_files\Mo2C\\orthorhombic\\ex_slab\\POSCAR").structure
#a choice of root slab that I have been using to test out the framework below



def generate_site_lst(slab, height=0.9):
    sf = AdsorbateSiteFinder(slab)
    #creates an AdsorbateSiteFinder object from pymatgen.analysis.adsorption 
    #to identify the possible sites for adsorbates    
    
    slab_corrected_surf = sf.assign_site_properties(slab, height=height)
    sf = AdsorbateSiteFinder(slab_corrected_surf)
    dict_ads_site = sf.find_adsorption_sites(distance=1.7, symm_reduce=False)
    #a dictionary with all possible adsorption sites on the surface
    
    site_lst = []
    
    for site_type in dict_ads_site.keys():
        i = 0
        if site_type != 'all':
#            if site_type != 'bridge':
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
        mol = Molecule(['C','O'],[np.array([0, 0, 0]),np.array([0, 0, 1.13])])
    else:
        mol = Molecule(['O','C'],[np.array([0, 0, 0]), np.array([0, 0, 1.13])])
    return mol    
    
def make_NH3():
    #creates an ammonia molecule to adsorb to the surface
    #only one orientation, because I've mostly only seen this configuration for NH3 on surfaces
    mol = Molecule(['N','H','H','H'], [np.array([0, 0, 0]),
                                                         np.array([-0.8248, -0.4762, 0.336]),
                                                         np.array([0.8248, -0.4762, 0.336]),
                                                         np.array([0, 0.9524, 0.336])])
    return mol

def make_H():
    #creates an hydrogen atom to adsorb to the surface
    mol = Molecule(['H'],[np.array([0,0,0])])
    return mol

def make_O():
    #creates an oxygen atom to adsorb to the surface
    mol = Molecule(['O'],[np.array([0,0,0])])
    return mol

def make_C():
    #creates an oxygen atom to adsorb to the surface
    mol = Molecule(['C'],[np.array([0,0,0])])
    return mol

def make_site_combos(slab,n_sites, height=0.9):
    #takes all of the possible adsorption site labels and makes all possible
    #combinations of a number n_sites of the sites with no repetition
    site_lst, dict_ads_site= generate_site_lst(slab, height=height)
    site_combos = []
    for combo in itertools.combinations(site_lst,n_sites):
        site_combos.append(combo)
    return site_combos

def determine_coverage(slab, coverage, ref_species=None, height=2.1):
    sf = AdsorbateSiteFinder(slab)
    surf_sites = sf.find_surface_sites_by_height(slab, height=height)
    if ref_species == None:
        n_surf_atoms = len(surf_sites)
    else:
        n_surf_atoms = 0
        for site in surf_sites:
            if site.species_string == ref_species:
                n_surf_atoms += 1
#    print(n_surf_atoms)
    n_sites_init = n_surf_atoms*coverage
    
    n_sites = np.round(n_sites_init)
    
    if n_sites != n_sites_init:
        actual_coverage = n_sites/n_surf_atoms
        print('Warning: the number of sites used does not exactly match the specified coverage, the actual coverage is %s' %actual_coverage)
    else:
        actual_coverage = coverage
    n_sites = int(n_sites)
    
    return n_sites, actual_coverage
    

def create_coord_combos(slab, coverage, ref_species=None, height=0.9, dist_reduce=2.1):
    
    n_sites, actual_coverage = determine_coverage(slab, coverage, ref_species=ref_species, height=height)
    site_combos = make_site_combos(slab,n_sites,height=height)
    #generates combinations of sites
    site_lst, dict_ads_site = generate_site_lst(slab, height=height)
    #generates individual sites
    idx_lst = []
    latt = slab.lattice.matrix
    #print(latt)
    ngh_latt_2d = [np.array([0,0,0]),np.array([1,0,0]),np.array([0,1,0]),
                       np.array([-1,0,0]), np.array([0,-1,0]),np.array([1,1,0]),
                       np.array([-1,1,0]), np.array([1,-1,0]),np.array([1,1,0])]
    ngh_latt_2d_metric = []
    for latt_vec in ngh_latt_2d:
        ngh_latt_2d_metric.append(latt_vec@latt)
    for idx in range(1, n_sites+1):
        neg_idx = idx * (-1)
        idx_lst.append(neg_idx)
    #a list of the indices for adsorbates used for calculating distance between adsorbates
    coord_combos = []
    for combo in site_combos:
        combo_dict = {}
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
            #stores a list of sites and appends hydrogen in those sites to perform distance measurements between the sites
        dirname = str(int(actual_coverage*100)) + 'ML/' + '_'.join(dirname_lst)
        combo_dict[dirname] = selected_sites
        dist_lst = []
        for idx_combo in itertools.combinations(idx_lst,2):
            dist = init_slab.get_distance(idx_combo[0], idx_combo[1])
            dist_lst.append(dist)
            #measures the distance between the sites
        if actual_coverage:
            if (n_sites == 1) or (dist_reduce < min(dist_lst)):
                coord_combos.append(combo_dict)
    return coord_combos

def save_site_combos(slab, adsorbate, path, coverage, height=0.9,
                     dist_reduce=2.1, symm_reduce=False, ref_species=None,no_bridge=False):
    coord_combos = create_coord_combos(slab, coverage, ref_species=ref_species, height=height, dist_reduce=dist_reduce)
    if symm_reduce:
        coord_combos = combo_symm_reduce(slab,coord_combos)
    if no_bridge:
        coord_combos = [name for name in coord_combos if 'bridge' not in list(name.keys())[0]]
    print(len(coord_combos))
    for combo in coord_combos:
        fin_slab = slab.copy()
        sites = list(combo.values())[0]
        dirname = list(combo.keys())[0]
        #sf = AdsorbateSiteFinder(fin_slab)
        for site in sites:
            sf = AdsorbateSiteFinder(fin_slab)
            fin_slab = sf.add_adsorbate(adsorbate,site,reorient=False)
        fin_slab = fin_slab.get_sorted_structure()
        #appends the specified adsorbate to the slab in the selected sites if the distance between the sites is more than a specified number of angstroms
            
        if not os.path.exists('%s\\%s'  %(path, dirname)):
            os.makedirs('%s\\%s' %(path, dirname))
        fin_slab.to('poscar','%s\\%s\\POSCAR'%(path, dirname))
        
        #stores the generated slabs with adsorbates in a local directory

def combo_symm_reduce(slab, coord_combos):
    surf_sg = SpacegroupAnalyzer(slab, 0.1)
    symm_ops = surf_sg.get_symmetry_operations()
    #get symmetry operations for the spacegroup of our slab
    test = []
    for coord in list(coord_combos[0].values())[0]:
        test.append(np.array([1,1,1]))
    unique_combos = [{'test_case':test}]
    #create a generic first combination of coordinates based on the number of 
    #coordinates in the combination submitted. necessary to initialize the search, 
    #but this set should never be matched and is removed at the end of the cycle
    u_combos_lst = []
    duplicate_lst = []
    for combo in coord_combos:
        coords = list(combo.values())[0]
        coords = [slab.lattice.get_fractional_coords(coord) for coord in coords]
        #converts coordinates to fractional form
        for u_combo in unique_combos:
            match = False
            u_coords = list(u_combo.values())[0]
            u_coords = [slab.lattice.get_fractional_coords(u_coord) for u_coord in u_coords]
            #converts unique coordinates in memory to fractional form
            for op in symm_ops:
                coord_match = []
                #iterates through symmetry operations
                for coord in coords:
                    incoord = False
                    if mg.util.coord.in_coord_list_pbc(u_coords, op.operate(coord), atol=1e-6):
                            #tests if the coordinate is in the combination of unique coordinates
                            incoord = True
                    coord_match.append(incoord)
                if all(coord_match):
                    #if all of the coordinates in a combination are already in the memory, 
                    #breaks the loop so that the repeated combination is not 
                    #added to the set of the unique coordinate combinations
                    #print('suggests %s is symmetrically equivalent to %s' %(list(combo.keys())[0], list(u_combo.keys()))[0])
#                    u_combos_lst.append(list(u_combo.keys())[0])
#                    duplicate_lst.append(list(combo.keys())[0])
                    match = True
                    break
            if match:
                break
        if not match:
            unique_combos.append(combo)
    unique_combos.pop(0)
#    combo_pairs = {'combo_pairs1':u_combos_lst, 'combo_pairs2':duplicate_lst}
#    pair_df = pd.DataFrame.from_dict(combo_pairs)
#    pair_df.to_csv(path_or_buf=os.getcwd()+'\\75ML_eq_pairs.csv')
    
                
    return unique_combos #, pair_df

def pair_compare(row, df):
    config1 = row['combo_pairs1']
    config2 = row['combo_pairs2']
    e1 = float(df['E'].get((df.site == './adsorbates/CO/clean/75ML_Cdown/' + config1)))
    e2 = float(df['E'].get((df.site == './adsorbates/CO/clean/75ML_Cdown/' + config2)))
    if abs(e1-e2) < 0.1:
        close_test = True
    else:
        close_test = False
    return e2
#    row['close'] = close_test

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
            slab = Poscar.from_file("C:\\Users\mjankous\Documents\VASP_files\Mo2C\\orthorhombic\\ex_slab\\POSCAR").structure
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
            slab = Poscar.from_file("C:\\Users\mjankous\Documents\VASP_files\Mo2C\\orthorhombic\\ex_slab\\POSCAR").structure
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
            slab = Poscar.from_file("C:\\Users\mjankous\Documents\VASP_files\Mo2C\\orthorhombic\\ex_slab\\POSCAR").structure
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

#def struc_to_graph(slab):
#    G = nx.Graph()
#    i = 0
#    for atom in slab:
#        i += 1
#        G.add_node(i, coords=atom.coords, species=atom.species_string)
#    G = NNedges(G)
#    return G
#
#def NNedges(G):
#    for nodei in G.nodes.data():
#            NN_dist = 100
#            NN_lst = []
#            for nodej in G.nodes.data():
#                if nodei != nodej:
#                    dist = distance.euclidean(nodei[1]['coords'],nodej[1]['coords'])
#                    if abs(dist - NN_dist) < 0.5:
#                        NN_lst.append(nodej[0])
#                    elif dist < NN_dist:
#                        NN_lst = [nodej[0]]
#                        NN_dist = dist
#            for node in NN_lst:
#                G.add_edge(nodei[0], node)
#    return G
#
#def periodic_isomorph(slab, ads_set1, ads_set2):
#    G = nx.Graph()
#    superslab = slab.copy()
#    superslab.make_supercell([2,2,1])
#    sf = AdsorbateSiteFinder(superslab)
#    surf_sites = sf.find_surface_sites_by_height(superslab,height=0.9)
#    i = 0
#    for site in superslab:
#        if site in surf_sites:
#            i += 1
#            G.add_node(i, coords=site.coords, species=site.species_string)
#    G1 = G.copy()
#    j=0
#    for adsorbate in ads_set1:
#        j += 1
#        G1.add_node('ads' + str(j), coords=adsorbate)
#    G2 = G.copy()
#    j=0
#    for adsorbate in ads_set2:
#        j += 1
#        G2.add_node('ads' + str(j), coords=adsorbate)
#    G1 = NNedges(G1)
#    G2 = NNedges(G2)
#    iso = nx.algorithms.isomorphism.is_isomorphic(G1,G2)
#    return iso, G1, G2
#    
#        
#def test_isomorph(slab1, slab2):
#    G1 = struc_to_graph(slab1)
#    G2 = struc_to_graph(slab2)
#    iso = nx.algorithms.isomorphism.is_isomorphic(G1,G2)
#    return iso
#

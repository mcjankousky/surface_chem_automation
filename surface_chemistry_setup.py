#!/usr/bin/env python

import argparse,move_files_in
import sitefind_tool as sft
import pandas as pd
from pymatgen.io.vasp.inputs import Poscar
import os

parser = argparse.ArgumentParser(description='This routine generates all possible configuration of a given adsorbate on the most stable surfaces of a material, given a minimum \
                                              distance allowed between molecules and a specified surface coverage. The defaults are a distance of 1.7 angstrom and a coverage of 0.5')
parser.add_argument('which_adsorbate', help='Which adsorbate molecule or atom is being put on the surface. Options are CO,NH3,H')
parser.add_argument('--logical_carbon_down', dest='logical_carbon_down',type=bool,default=True,
                    help='If adsorbate is CO, logical giving whether carbon atom is closest to surface. Default is True, i.e. carbon is closest to surface')
parser.add_argument('--run_directory', dest='run_dir',default='./',help='Location in which to run this script.')
parser.add_argument('--number_surf', dest='number_surf',type=int,default=4,help='Number of surfaces for which to perform surface chemistry calculations')
parser.add_argument('--min_ads_spacing', dest='min_ads_spacing',type=float,default=1.7,help='Minimum spacing between adsorbates on a surface, in angstrom. Default is 1.7 angstrom.')
parser.add_argument('--coverage', dest='coverage',type=float,default=0.5,help='Coverage of molecules on actives sites on the surfaces. Default is 0.5')
parser.add_argument('--phase_evaluated', dest='phase_evaluated', type=bool, default=False, help='Labels whether or not to use the stability csv that has been sorted based on average surface energy in a given chemical potential')
parser.add_argument('--ref_species', dest='ref_species',type=str, default=None, help='The species at the surface to which adsorbate coverage is compared')
args = parser.parse_args()

run_dir=args.run_dir
number_surf=args.number_surf
phase_evaluated = args.phase_evaluated
# First read in the surface_stability csv, and take top number_surf surfaces on which to do surface chemistry calculations
if not phase_evaluated:
    E_df=pd.read_csv(run_dir + 'surface_stability.csv',sep=',')
else:
    E_df=pd.read_csv(run_dir + 'surface_stability_ave_E_sorted.csv',sep=',')
surf_names=E_df['calculation']
print(surf_names.tolist())
try:
    surf_names=surf_names[0:number_surf]
except:  # DVF comment: I should be more specific here, but this is meant to handle if there is no data in the surface_stability.csv file. 
    surf_names=surf_names  

if args.which_adsorbate=='CO':
    mol=sft.make_CO(args.logical_carbon_down)
elif args.which_adsorbate=='NH3':
    mol=sft.make_NH3()
else:
    mol=sft.make_H()

symm_reduce=False  #DVF: this can be added as an input flag if we ever figure this out


for surf in surf_names.tolist():
    surf_path=run_dir + 'surface_stability/'+surf
    print(surf_path)
    slab = Poscar.from_file(surf_path+'/POSCAR').structure
    site_list,dict_ads_site=sft.generate_site_lst(slab)
    surf_chem_path=run_dir+'/surface_chemistry/'+ surf + '/' + args.which_adsorbate + '/'
    sft.save_site_combos(slab, mol, surf_chem_path, args.coverage, args.min_ads_spacing, symm_reduce, args.ref_species)   # DVF comment: I am converting coverage to an int here
                                                                                                                 # so n_sites can still be passed here while you transition this over.                                                                                     
                                                                                                                 
    if not os.path.exists(run_dir + '/surface_chemistry/' + surf + '/reference_surface'):
        os.makedirs(run_dir + '/surface_chemistry/' + surf + '/reference_surface')
    slab.to("poscar",run_dir+'/surface_chemistry/' + surf + '/reference_surface/POSCAR')

move_files_in.move_files_in(run_dir,'surface_chemistry')                                                                                


#!/usr/bin/env python

import argparse,E_conv_collect#,phase_diagram DVF comment: phase diagram routine should be written
import pandas as pd

parser = argparse.ArgumentParser(description='This routine collects total energies for a set of surfaces computed by VASP, and outputs a csv with these energies, and a phase diagram (in development).')
parser.add_argument('--run_directory', default='./',help='Location in which to run this script.')
parser.add_argument('--phase_diagram', type=bool,default=False,help='Whether or not to generate a phase-stability plot. Default is false.')
parser.add_argument('--number_surf', type=int,default=4,help='Number of surfaces for which to generate phase diagram plot. Default is 4.')
parser.add_argument('--chemical_potential_range', type=float, default=[0,10], nargs=2,metavar=('chem_pot_low', 'chem_pot_high'),
                    help='Range of chemical potentials for phase diagram plots. Defaults to 0-10.')
args = parser.parse_args()

run_dir=args.run_directory
# Extract energies from all surface calculations, put them in pandas 
# data frame and write as csv. 
surf_dir=run_dir+'surface_stability'
E_df = E_conv_collect.E_conv_collect(surf_dir,'surface_stability')
# If specified, generate phase diagram
if args.phase_diagram:
    phase_diagram.create_phase_diagram_surf(surf_dir,number_surf,E_df,args.chemical_potential_range)



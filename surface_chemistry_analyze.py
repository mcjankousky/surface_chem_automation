#!/usr/bin/env python

import argparse,E_conv_collect,phase_diagram #DVF comment: phase diagram routine should be written
import pandas as pd

parser = argparse.ArgumentParser(description='This routine collects total energies for all configurations of an adsorbates on the stable surfaces of a bulk material computed by VASP, \
                                              and outputs a csv with these energies, or a phase diagram (in development).')
parser.add_argument('--run_directory', dest='run_dir',default='./',help='Location in which to run this script.')
parser.add_argument('--phase_diagram', dest='phase_diagram',type=bool,default=False,help='Whether or not to generate a phase-stability plot. Default is false.')
parser.add_argument('--number_coverage', dest='num_coverage',type=float,default=4,help='Number of coverages for which to generate phase diagram plot. Default is 4.')
parser.add_argument('--chemical_potential_range', type=float, default=[0,10], nargs=2,metavar=('chem_pot_low', 'chem_pot_high'),
                    help='Range of chemical potentials for phase diagram plots. Defaults to 0-10.')
parser.add_argument('--pressure_range', type=float, default=[-4,-0.01], nargs=2,metavar=('pressure_low', 'pressure_high'),
                    help='Range of natural log of species partial pressure for phase diagram plots. Defaults to -4-(-0.01).')
parser.add_argument('--temperature_range', type=float, default=[1,500], nargs=2,metavar=('temp_low', 'temp_high'),
                    help='Range of temperatures for phase diagram plots in Kelvin. Defaults to 1-500.')
args = parser.parse_args()

run_dir=args.run_dir
# Extract energies from all surface calculations, put them in pandas 
# data frame and write as csv. 
chem_dir=run_dir#+'surface_chemistry'
E_df = E_conv_collect.E_conv_collect(chem_dir,'surface_chemistry',False)
# If specified, generate phase diagram
if args.phase_diagram:
    phase_diagram.create_phase_diagram_chem(chem_dir,args.num_coverage,E_df,args.pressure_range,args.temperature_range) 



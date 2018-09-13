#!/usr/bin/env python

import argparse,slab_factory,move_files_in
from pymatgen.io.vasp.inputs import Poscar

parser = argparse.ArgumentParser(description='This routine generates all possible surface terminations for a given bulk material (described in a POSCAR) and maximum miller index, and \
                                              sets up directories in which VASP calculations can be run to determine relative surface stabilities.')
parser.add_argument('max_miller_ind',type=int,help='maximum miller index for surface terminations')
parser.add_argument('slab_thickness',type=float,
                    help='thickness of slab to generate in angstrom. Slabs within +-3 angstrom of this thickness will be considered')  # DVF: could add tolerance as an optional variable
parser.add_argument('--surface_supercell', dest='surf_supercell',type=int, default=[1,1], nargs=2,metavar=('nss1', 'nss2'),help='surface supercell dimensions in each surface lattice direction. Defaults to 1.')
parser.add_argument('--run_directory', dest='run_dir',default='./',help='Location in which to run this script.')
args = parser.parse_args()

run_dir=args.run_dir

struc = Poscar.from_file(run_dir + 'POSCAR').structure

# DVF comment: the slab_factory routine uses the range function with slab_thickness, hence
# the integer conversion here. I doubt this is a big deal. We may want to consider using
# in_unit_planes=True so that we do things in terms of numbers of layers, which I think
# may be more intuitive to users. 
slab_list = slab_factory.get_all_slabs(struc,args.max_miller_ind,int(args.slab_thickness),args.surf_supercell,run_dir)

move_files_in.move_files_in(run_dir,'surface_stability')


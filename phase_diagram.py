# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 11:53:33 2018

@author: mjankous
"""

import os
import pandas as pd
import ast
import numpy as np
import matplotlib.pyplot as plt

def create_phase_diagram_surf(surf_dir,number_surf,E_df,chemical_potential_species,
                              chemical_potential_range, plot_title='Surface Stability Phase Diagram'):
    #E_df['species'] = E_df.apply(lambda row: species_convert(row), axis=1)
    for row in E_df.iterrows():
        if '/bulk' in row[1]['calculation']:
            bulk_row = row[1]
            break
    bulk_E = bulk_row['E']
    bulk_species = ast.literal_eval(bulk_row['species'])
    E_df['surf_E'] = E_df.apply(lambda row: calculate_surf_E(row, bulk_E, 
                                                            bulk_species, chemical_potential_species, 
                                                            chemical_potential_range), axis = 1)
    E_df['average_surf_E'] = E_df.apply(lambda row: calculate_ave_surf_E(row), axis = 1)
    E_df.sort_values(by='average_surf_E', ascending=True, inplace=True)
    mu_ref = np.linspace(chemical_potential_range[0],chemical_potential_range[1])
    plot_generator = phase_diagram_plotmake()
    plot_generator.make_surf_stability_phase(E_df,mu_ref,number_surf, plot_title, chemical_potential_species, surf_dir)
    E_df.to_csv(surf_dir + '/surface_stability_ave_E_sorted.csv')
    #return E_df

def calculate_surf_E(row, bulk_E, bulk_species, chemical_potential_species, chemical_potential_range):
    row_species = ast.literal_eval(row['species'])
    if len(bulk_species) == 1:
        Nbulk = 0
        Nslab = 0
        for bulk_value in bulk_species.values():
            Nbulk += bulk_value
        for slab_value in row_species.values():
            Nslab += slab_value
        surf_E = 0.5*(row['E'] - (Nslab/Nbulk)*bulk_E)
    elif len(bulk_species) == 2:
        try:
            num_ref_species_slab = row_species[chemical_potential_species]
        except KeyError:
            num_ref_species_slab = 0
        num_ref_species_bulk = bulk_species[chemical_potential_species]
        for bulk_key in bulk_species.keys():
            if bulk_key != chemical_potential_species:
                other_species = bulk_key
        try:
            num_other_species_slab = row_species[other_species]
        except KeyError:
            num_other_species_slab = 0
        num_other_species_bulk = bulk_species[other_species]
        mu_ref = np.linspace(chemical_potential_range[0], chemical_potential_range[1])
        mu_other = bulk_E/num_other_species_bulk - (num_ref_species_bulk/num_other_species_bulk)*mu_ref
        surf_E = (row['E'] - num_other_species_slab*mu_other - num_ref_species_slab*mu_ref)/float(row['area'])
        surf_E = list(surf_E)        
    return surf_E

def calculate_ave_surf_E(row):
    average_surf_E = np.average(row['surf_E'])
    return average_surf_E

def species_convert(row):
    species_dict = ast.literal_eval(row['species'])
    return species_dict

class phase_diagram_plotmake(object):
    
    def fig_size(self, fig_size=None):
        fig_size = (fig_size or (9.2,7))
        return fig_size
    
    @property
    def tableau20(self):
    
        tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
        
                     (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
        
                     (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
        
                     (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
        
                     (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    
        # Rescale to values between 0 and 1 
        
        for i in range(len(tableau20)):  
        
            r, g, b = tableau20[i]  
        
            tableau20[i] = (r / 255., g / 255., b / 255.)
        
        
        
        # Use with plt.plot(…, color=tableau[0],…)
        return tableau20
    
    def make_surf_stability_phase(self,E_df,mu_ref,number_surf, plot_title, chemical_potential_species, surf_dir):
        fig, ax = plt.subplots(figsize=(23,15))
        plt.xticks(fontsize=36)
        N = 0
        for row in E_df.iterrows():
            calc = row[1]['calculation']
            surf_E = row[1]['surf_E']
            if N < 10:
                ax.plot(mu_ref,surf_E, linestyle='-', color=self.tableau20[2*N], linewidth=4, label = calc[2:])
            elif N < 20:
                ax.plot(mu_ref,surf_E, linestyle='-',color=self.tableau20[2*(N-10)+1], linewidth=4, label=calc[2:])
            N += 1
            if N > number_surf:
                break
        ax.tick_params(axis = 'both', which = 'major', labelsize=40)
        ax.set_xlabel("$\mu_%s$ (eV)" %chemical_potential_species,size=44)
        ax.set_ylabel(r"Surface Energy (eV/$\mathring{A}^2$)",size=44)
        ax.set_title(plot_title,size=60)
        plt.legend(fontsize=36)
        #plt.show()
        fig.savefig(surf_dir + '/phase_diagram.png')

def specify_species(row):
    if '/H/' in row['site']:
        return "H"
    elif '/NH3/' in row['site']:
        return "NH3"
    elif '/CO/' in row['site']:
        return "CO"
    elif '/O_coverage/' in row['site']:
        return "O"
    elif '/O/' in row['site']:
        return "O"
    else:
        return None
    
def specify_ref_surf(row):
    if '/clean/' in row['site']:
        return 'clean'
    elif '/5ML_O/' in row['site']:
        return '5ML_O'
    elif 'gas_phase_ref' in row['site']:
        return 'gas_phase_ref'
    else:
        return 'clean'

def specify_coverage(row):
    site_data = row['site'].split('/')
    for item in site_data:
        if "ML" in item:
            if item != '5ML_O':
                coverage = item
                break
        else:
            coverage = None
    return coverage
        
def specify_config(row):
    site_data = row['site'].split('/')
    return site_data[-1]

def make_false_0(row):
    if row['conv'] == False:
        return 0
    else:
        return row['E']

def conv_E_to_float(row):
    return float(row['E'])

ev_per_j = 6.241506e18

atom_per_mol = 6.022e23



thermo_constants = {}

thermo_constants['H2'] = {}
thermo_constants['H2']['A'] = 33.066178
thermo_constants['H2']['B'] = -11.363417
thermo_constants['H2']['C'] = 11.432816
thermo_constants['H2']['D'] = -2.772874
thermo_constants['H2']['E'] = -0.158558
thermo_constants['H2']['F'] = -9.980797
thermo_constants['H2']['G'] = 172.707974
thermo_constants['H2']['H'] = 0
thermo_constants['H2']['S'] = 130.68
thermo_constants['H2']['E_dft'] = -6.76 #eV per atom
thermo_constants['H2']['ref_to_ads'] = 0.5

thermo_constants['O2'] = {}
thermo_constants['O2']['A'] = 31.32234
thermo_constants['O2']['B'] = -20.23531
thermo_constants['O2']['C'] = 57.86644
thermo_constants['O2']['D'] = -36.50624
thermo_constants['O2']['E'] = -0.007374
thermo_constants['O2']['F'] = -8.9033471
thermo_constants['O2']['G'] = 246.7945
thermo_constants['O2']['H'] = 0
thermo_constants['O2']['S'] = 205.15
thermo_constants['O2']['E_dft'] = -9.85 #eV per atom
thermo_constants['O2']['ref_to_ads'] = 0.5

thermo_constants['NH3'] = {}
thermo_constants['NH3']['A'] = 19.99563
thermo_constants['NH3']['B'] = 49.77119
thermo_constants['NH3']['C'] = -15.37599
thermo_constants['NH3']['D'] = 1.921168
thermo_constants['NH3']['E'] = 0.189174
thermo_constants['NH3']['F'] = -53.30667
thermo_constants['NH3']['G'] = 203.8591
thermo_constants['NH3']['H'] = -45.89806
thermo_constants['NH3']['S'] = 192.77
thermo_constants['NH3']['E_dft'] = -19.5 #eV per atom
thermo_constants['NH3']['ref_to_ads'] = 1

thermo_constants['CO'] = {}
thermo_constants['CO']['A'] = 25.56759
thermo_constants['CO']['B'] = 6.096130
thermo_constants['CO']['C'] = 4.054656
thermo_constants['CO']['D'] = -2.671301
thermo_constants['CO']['E'] = 0.131021
thermo_constants['CO']['F'] = -118.0089
thermo_constants['CO']['G'] = 227.3665
thermo_constants['CO']['H'] = -110.5271
thermo_constants['CO']['S'] = 197.66
thermo_constants['CO']['E_dft'] = -14.8 #eV per atom
thermo_constants['CO']['ref_to_ads'] = 1

kb = 1.380648*(10**(-23))


def create_phase_diagram_chem(chem_dir,coverage_values,E_df,ads_species,ref_surf,ref_species,temperature_range,pressure_range):
    if ref_surf == 'clean':
        for row in E_df.iterrows():
            if '/clean_surface' in row[1]['site']:
                clean_row = row[1]
                break
        clean_E = float(clean_row['E'])
    elif ref_surf == '5ML_O':
        for row in E_df.iterrows():
            if '/O_coverage/5ML/hollow_combos/12' in row[1]['site']:
                clean_row = row[1]
                break
        clean_E = float(clean_row['E'])
    E_df['species'] = E_df.apply(lambda row: specify_species(row), axis=1)
    E_df['ref_surf'] = E_df.apply(lambda row: specify_ref_surf(row), axis=1)
    E_df['coverage'] = E_df.apply(lambda row: specify_coverage(row), axis=1)
    E_df['configuration'] = E_df.apply(lambda row: specify_config(row), axis=1)
    E_df['E'] = E_df.apply(lambda row: make_false_0(row), axis=1)
    E_df['E'] = E_df.apply(lambda row: conv_E_to_float(row), axis=1)
    e_dict = {}
    for coverage in coverage_values.keys():
        print(coverage, ref_surf, ads_species)
        e_dict[coverage] = min(E_df.E.get((E_df['species'] == ads_species) & (E_df['ref_surf'] == ref_surf) & (E_df['coverage'] == coverage)).tolist())
    
    yseed = np.linspace(temperature_range[0], temperature_range[1], 200)
    xseed = np.linspace(pressure_range[0], pressure_range[1], 200)
    karray = []
#    kmax = max(coverage_values.values())
#    print(kmax)
    for p in xseed:
        kx = []
        for T in yseed:
            dH_gas = 1000*(thermo_constants[ref_species]['A']*(T/1000) + thermo_constants[ref_species]['B']*((T/1000)**2)/2 + thermo_constants[ref_species]['C']*((T/1000)**3)/3 + thermo_constants[ref_species]['D']*((T/1000)**4)/4 - thermo_constants[ref_species]['E']/(T/1000) + thermo_constants[ref_species]['F'] - thermo_constants[ref_species]['H'])/(atom_per_mol)
            E_gas = thermo_constants[ref_species]['E_dft']/ev_per_j + dH_gas
            S_gas = (thermo_constants[ref_species]['A']*np.log(T/1000) + thermo_constants[ref_species]['B']*(T/1000) + thermo_constants[ref_species]['C']*((T/1000)**2)/2 + thermo_constants[ref_species]['D']*((T/1000)**3)/3 - thermo_constants[ref_species]['E']/(2*(T/1000)**2) + thermo_constants[ref_species]['G'])/(atom_per_mol)
            mu_gas = kb * T * p + E_gas - T * S_gas
            omega = {}
            for coverage in coverage_values.keys():
                omega[coverage] = e_dict[coverage]/ev_per_j - clean_E/ev_per_j - thermo_constants[ref_species]['ref_to_ads'] * coverage_values[coverage] * mu_gas
            if min(omega.values()) > 0:
                k = 0#kmax
            else:
                for coverage in coverage_values.keys():
                    if omega[coverage] == min(omega.values()):
                        k = coverage_values[coverage]#*(-1) + kmax
                        break
                    else:
#                        print(omega[coverage], min(omega.values()))
                        k = ""
#            print(T,p,omega)
            kx.append(k)
        karray.append(kx)
    
    karray = np.asarray(karray)

    karray = karray.transpose()
#    fig,ax = plt.subplots(figsize=(15,8))
    plt.figure(figsize=(15,8))
    #plt.imshow(ki, vmin=min(k_lst), vmax=max(k_lst), origin='lower',
    #           extent=[min(lnpp_lst), max(lnpp_lst), min(T_lst), max(T_lst)],
    #           aspect='auto', cmap='winter')
    plt.contourf(xseed, yseed, karray, 20, cmap='winter')
    plt.ylabel('T (K)', fontsize=36)
    plt.yticks(fontsize=28)
    plt.xlabel(r'$ln(\frac{p_{%s}}{p_o})$' %ref_species, fontsize=50)
    plt.xticks(fontsize=28)
    plt.title(r'%s Coverage on Surface' %ads_species, fontsize=36)
    plt.colorbar()
    plt.show()
    plt.savefig(chem_dir + '/' + ads_species + '_coverage_diagram.png')
    
#    return E_df
    
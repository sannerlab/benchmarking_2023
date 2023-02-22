#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 15:03:30 2022
@author: sshanker

This code is used to make Venn diagrams for accurately predicted complexes
from all methods
"""

import pandas as pd
import numpy as np
import itertools

Fnat_cut_off = 0.8
cross_analysis_between = ['AF2_multi_5_Fnat', 'ADCP_5_Fnat', 'AF2_mono_5_Fnat','OF_1_Fnat']
cross_analysis_between = ['AF2_multi_1_Fnat', 'ADCP_1_Fnat', 'AF2_mono_1_Fnat','OF_1_Fnat']
# cross_analysis_between = ['AF2_multi_5_Fnat', 'ADCP_30_Fnat']



def generate_combinations(length):
    '''
    To generate all combinations of given methods
    '''
    out = []
    
    for i in range(1,length+1):
        for itr in itertools.combinations(range(length), i):
            out.append(list(itr))
    
    return out


def calculate_venn_values():
    fnat_file = 'all_fnat.dat'
    fnat_data = pd.read_csv(fnat_file,sep = ",\s+")

    all_index_above_cutoff=[]
    for i in fnat_data.columns:  
        if i == 'pdb':
            continue
        all_index_above_cutoff.append(np.where(fnat_data[i] >= Fnat_cut_off)[0])
    
    columns = list(fnat_data.columns)

    all_indexes =[]
    for method_ in cross_analysis_between:
        all_indexes.append(all_index_above_cutoff[columns.index(method_) -1])


    max_v = 0
    for  i in all_indexes:
        max_v = max(max_v, max(i))
    
    pos_array_of_methods = [[] for i in range(max_v+1)]
    for i, ind in enumerate(all_indexes):
        for val in ind:
           pos_array_of_methods[val].append(i)         
    
    for m in generate_combinations(len(cross_analysis_between)):      
        print_seq = ''
        for ml in m:
            print_seq = (print_seq + 
                         cross_analysis_between[ml].replace('Fnat','' ).replace("_",'') +
                         ' & ')
        
        print(print_seq[:-2],' = ', pos_array_of_methods.count(m))
        

calculate_venn_values()




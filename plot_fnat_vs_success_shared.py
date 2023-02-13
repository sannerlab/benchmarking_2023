#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 12:09:24 2022
@author: sshanker
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


colors=['#1f37f9', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'] # manual colors for plot

# colors arranged to use per method
use_colors = [colors[0], colors[0], colors[0],              # for AF2 multi
              colors[1], colors[1],                         # for AF2 mono
              colors[3],                                    # for OF
              colors[2], colors[2] ,colors[2] ,colors[2]]   # for ADCP

# line style arranged to use per method
use_style =["", "--" , ".-",                                # for AF2 multi
            "", "--",                                       # for AF2 mono
            "",                                             # for OF
            "","--", ".-",":",]                             # for ADCP

fnat_file = 'all_fnat.dat' # fnat data file for all pdb
data_ = pd.read_csv(fnat_file,sep = ",\s+")

histograms = np.array(range(101))/100
value_array = []

for h in histograms:
    vals = np.array(np.sum(data_[
    ['AF2_multi_1_Fnat', 'AF2_multi_5_Fnat', 'AF2_multi_25_Fnat',
     'AF2_mono_1_Fnat', 'AF2_mono_5_Fnat',
     'OF_1_Fnat',
     'ADCP_1_Fnat','ADCP_5_Fnat', 'ADCP_30_Fnat', 'ADCP_all_Fnat']
    ]>=(h),0)) #to match excel sheet

    value_array.append(vals)
    
value_array = np.array(value_array)
value_array = 100*value_array/value_array[0,:]

plt.subplots(1,1, figsize =(16,11))
tick_size = 20
legend_size=13.5


# Only ADCP plots
plt.subplot(2,2,1)
for i in range(4):
    plt.plot(histograms, value_array[:,6+i], 
             use_style[6+i],
             color=use_colors[6+i],
             markersize =7)
    
plt.legend(['ADCP 1', 'ADCP 5', 'ADCP 30', 'ADCP All',], fontsize=legend_size)
plt.plot([0.8,0.8],[0,100],'k--')
plt.xlim([1,0])
plt.ylim([-1,101])

plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
plt.ylabel('% Complexes',fontsize=25)

# ADCP + OF plots
plt.subplot(2,2,2)
for i in range(5):
    plt.plot(histograms,value_array[:,5+i], 
             use_style[5+i],
             color=use_colors[5+i],
             markersize =7)
    
plt.legend(['OF',
            'ADCP 1', 'ADCP 5', 'ADCP 30', 'ADCP All',],
           fontsize=legend_size,loc=4)

plt.plot([0.8,0.8],[0,100],'k--')
plt.xlim([1,0])
plt.ylim([-1,101])
# plt.xlabel(r'Fraction of Native Contacts',fontsize=18)
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)

#ADCP + AF2mono plots
plt.subplot(2,2,3)
for i in range(7):
    if i==2:
        continue # ignore OF data
    plt.plot(histograms,value_array[:,3+i], 
             use_style[3+i],
             color=use_colors[3+i],
             markersize =7)
    
plt.legend(['AF2mono 1', 'AF2mono 5',
            'ADCP 1', 'ADCP 5', 'ADCP 30', 'ADCP All'],
           fontsize=legend_size)

plt.plot([0.8,0.8],[0,100],'k--')
plt.xlim([1,0])
plt.ylim([-1,101])
plt.xlabel(r'Fraction of Native Contacts',fontsize=25)
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
plt.ylabel('% Complexes',fontsize=25)

# ADCP + AF2multi plots
plt.subplot(2,2,4)
for i in range(10):
    if [3,4,5].count(i)>0:
        continue # ignore OF and AF2 mono data
    plt.plot(histograms,value_array[:,i], 
             use_style[i],
             color=use_colors[i],
             markersize =7)
    
plt.legend(['AF2multi 1', 'AF2multi 5', 'AF2multi 25',
            'ADCP 1', 'ADCP 5', 'ADCP 30', 'ADCP All',],
           fontsize=legend_size)

plt.plot([0.8,0.8],[0,100],'k--')
plt.xlim([1,0])
plt.ylim([-1,101])
plt.xlabel(r'Fraction of Native Contacts',fontsize=25)
plt.xticks(fontsize=tick_size)
plt.yticks(fontsize=tick_size)
plt.tight_layout(w_pad=1)


# Manually adjusted subplot lables
plt.text(2.13,200,'A', fontsize =45,weight="bold")
plt.text(0.98,200,'B', fontsize =45,weight="bold")
plt.text(2.13,85,'C', fontsize =45,weight="bold")
plt.text(0.98,85,'D', fontsize =45,weight="bold")
# plt.ylabel('% Complexes',fontsize=18)

plt.savefig('Fnat2.png',dpi=300)


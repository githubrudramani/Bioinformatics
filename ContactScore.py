#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 10:46:34 2020

@author: rudramani

It is a simple model to find the contact score between two proteins in 
contact with each other. It requires residues pairs between proteins in 
contact.

"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.spatial import distance
from sklearn.preprocessing import StandardScaler




polar = ['SER', 'THR', 'TYR', 'HSD', 'CYS', 'ASN', 'GLN']
nonpolar = ['ALA','VAL','PHE','ILE','LEU','PRO','MET','TRP','GLY']
chargedPlus = ['LYS', 'ARG']
chargedMinus = ['ASP', 'GLU']

Weights = {'polar_polar': -1, 'polar_nonpolar':0, 'polar_chargedPlus':-2, \
        'polar_chargedMinus':2,'nonpolar_nonpolar':3, 'nonpolar_chargedPlus':0, \
        'nonpolar_chargedMinus':0, 'chargedPlus_chargedPlus':-4,\
        'chargedPlus_chargedMinus':4, 'chargedMinus_chargedMinus':-4, }


##### Two kind of class object to find restype
'''
class ResType:
    def f(self,a):
        if a in polar:
            self.type = "polar"
        elif a in nonpolar:
            self.type = "nonpolar"
        elif a in chargedPlus:
            self.type = "chargedPlus"
        else:
            self.type = "chargedMinus"

x = ResType()
x.f('GLU')
x.type
'''

class ResType2:
    def __init__(self,a):
        if a in polar:
            self.type = "polar"
        elif a in nonpolar:
            self.type = "nonpolar"
        elif a in chargedPlus:
            self.type = "chargedPlus"
        else:
            self.type = "chargedMinus"

        
#y = ResType2('HSD')
#y.type
            
###  Function to find Restype using restype class
def kind(resname):
    y = ResType2(resname)
    return y.type

##### Function to find contact score from residue pair
def ContactScore(Dataframe):
    score = 0
    for index, row  in fin.iterrows():
        pair = row.values
        pair1 = kind(pair[0])+'_'+kind(pair[1])
        pair2 = kind(pair[1])+'_'+kind(pair[0])
        if pair1 in Weights:
            score += Weights[pair1]
        elif pair2 in Weights:
            score += Weights[pair2]
    return score
    

### reading the residue pair within 3.5 of each other created by tcl script Pair.tcl
data = 'Pair_list.txt'
fin = pd.read_csv(data,delimiter='\s+', header=None)
columns = fin.iloc[0].values
fin = fin[1:]

score = ContactScore(fin)

print('contact score of this pdb is:', score)
        
    
    


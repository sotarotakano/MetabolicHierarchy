#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.0'

'''The codes for computing ranks in the growth-rate of a given set of
carbon sources in a given models. We specifically focued on the CarveMe
prokaryotic models and computed their ranks, but the algorithm is compartible with 
xml.gz models with similar format of CarveMe and BiGG database'''

import routine_functions.CAFBA_mod as CF
from routine_functions.computing_growth_flux import *
from pathlib import Path
import os
from os.path import join
import pandas as pd
import numpy as np
import cobra
from multiprocessing import Pool
import gzip
import copy
import glob
import gc
from itertools import combinations


def load_xml_gz_model(filepath):
    model = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(gzip.open(filepath,'rt')))
    model.reassign_compartments()
    return model

def load_xml_model(filepath):
    model = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(filepath))
    model.reassign_compartments()
    return model

def rank_data_matrix(dataframe,min_threshold = 1e-05,non_zero = False):
    #from scipy.stats import rankdata
    if non_zero:
        dataframe = dataframe.loc[:,dataframe.iloc[-1,:] > 0]
    Carbon_rank = copy.copy(dataframe)
    #index_list = dataframe.index
    for N in range(0,dataframe.shape[0]):
        growth_data = copy.copy(dataframe.iloc[N,:]).to_numpy()
        growth_data[growth_data < 1e-06] = -1
        rank_data = np.zeros(len(growth_data))
        i, k = 1, 0
        while max(growth_data) > 1e-06:
            max_growth = max(growth_data)
            rank_data[growth_data == max_growth] = i
            k += len(growth_data[growth_data == max_growth])
            growth_data[growth_data == max_growth] = 1e-06
            if np.log10(max_growth/max(growth_data)) > min_threshold:
                i = i + k
                k = 0
        if min(growth_data) < 0:
            rank_data[growth_data < 0] = max(rank_data) + 1
        Carbon_rank.iloc[N,:] = rank_data
    return Carbon_rank


def check_growth(arguments):
    filepath = arguments[0]
    savefile = arguments[1]
    if os.path.splitext(filepath)[-1] == ".gz":
        model = load_xml_gz_model(filepath)
    elif os.path.splitext(filepath)[-1] == ".xml":
        model = load_xml_model(filepath)
    else:
        return
    Growth = computing_growth(Sugar_list, model, name = model.id, 
                              lb = -120.0, normalize_by_Catoms = True, ErrorMessage = False)
    dataline = [str(x) for x in Growth]
    with open(savefile, mode = "a") as f:
        f.write(",".join([model.id] + dataline) + "\n")
    del model
    del Growth
    gc.collect()

#################################
#### Main body of the script ####
#################################

# Parameters
homedir = Path(__file__).parent.absolute()
cobra.Configuration().solver = "cplex" #setting optimizer
Nthreads = 6 #The number of threads

'''Here, the minimum number of models are specified. The program searching for the set of Csources
supporting the growth of > "min_models" '''
min_models = 80

# The list of sugars used here is same as Takano et al., 2023
Sugar_list = ['EX_glc__D_e', 'EX_fru_e', 'EX_man_e','EX_fuc__L_e', 'EX_melib_e','EX_gal_e', 'EX_rib__D_e']

# Loading models from the directory for CarveMe models. It's originally from Machado et al., 2018
model_dir = join(homedir,"carveme_mastermodels")
#taxon_info = pd.read_csv(join(home_dir,"auxiliary","Model_taxon.csv"), index_col = 0)
carveme_models = [x for x in glob.glob(join(model_dir,"*","*","*"+".xml.gz")) if "cafba" not in x]

print("[INPUTS] We found %i carveme models."%len(carveme_models))
print("[RUNNING] Checking the growth capacity...")

'''First, we checked the growth capability in the given all models on a given set of substrates.'''
growth_result = join(model_dir,"growth_carvememodels.csv")
if not os.path.exists(growth_result):
    with open(growth_result, mode = "w") as f:
        f.write(",".join(["model"] + Sugar_list) + "\n")
    
    for i in range(int(len(carveme_models)/(Nthreads*4))):
        with Pool(processes=Nthreads) as pool:
            pool.map(check_growth, [(x,growth_result) for x 
                                    in carveme_models[i*Nthreads*4:min(len(carveme_models)-1,int((i+1)*Nthreads*4))]])

carveme_growth = pd.read_csv(growth_result,index_col = 0)

'''Then we searched for the set of substrates satisfying the constraint 
(i.e., more than "min_models" of models can metabolize all of them).'''
max_grow = 1
for i in [int(x) for x in np.linspace(len(Sugar_list)-1,0,len(Sugar_list))]:
    for s in combinations(Sugar_list,i):
        grow = len(carveme_growth.loc[carveme_growth[list(s)].min(axis=1) > 1e-06,:])
        if grow > max_grow:
            max_grow = grow
            selected_sugars = s
    if max_grow >= min_models:
        break
    #print("%i models can grow on %s"%(grow,s))

print("The following combination of sugars supports the growth of %i models"%max_grow)
print(selected_sugars)
selected_growth = carveme_growth.loc[carveme_growth[list(selected_sugars)].min(axis=1) > 1e-06,selected_sugars]


'''Then we generates CAFBA models that can grow on a selected set of substrates.'''
print("[RUNNING] Creating CAFBA models for the selected CarveMe models...")
#convert to CAFBA and simulate
cm_valid_models = [x for x in carveme_models 
                   if len([y for y in selected_growth.index if y.replace("_xml","") in x]) > 0
                   and "cafba" not in x]
print("%i models can grow on those substrates and covert to CAFBA format..."%len(cm_valid_models))
carveme_models_CF = []
for cm in cm_valid_models:
    cafba_file_name = join("/".join(cm.split("/")[:-1]),cm.split("/")[-1].split(".")[0] + "_cafba.xml.gz")
    if not os.path.exists(cafba_file_name):
        CF.main(cm,output_file=cafba_file_name)
    if os.path.exists(cafba_file_name):
        carveme_models_CF.append(cafba_file_name)

# Load carveme cafba model
growth_result_cafba = join(model_dir,"growth_carvememodels_CAFBA.csv")
print("[RUNNING] Computing the growth-rate ranks for the selected models...")

if not os.path.exists(growth_result_cafba):
    with open(growth_result_cafba, mode = "w") as f:
        f.write(",".join(["model"] + Sugar_list) + "\n")

    # Check the growth rate of CAFBA models
    with Pool(processes=Nthreads) as pool:
        pool.map(check_growth, [[x, growth_result_cafba] for x in carveme_models_CF])

carveme_growth_CF = pd.read_csv(growth_result_cafba,index_col = 0)
carveme_growth_CF = carveme_growth_CF.loc[selected_growth.index,selected_sugars]

# Then computing the growth-ranks from the growth-rate data
carveme_rank_CF = rank_data_matrix(carveme_growth_CF)

#print("[RUNNING] Adding taxonomic informations...")
# Adding the taxon info. We used NCBI format here.
#for i in carveme_rank_CF.index:
#    carveme_rank_CF.loc[i,"taxon"] = taxon_info.loc[i,"taxon"]
#    carveme_rank_CF.loc[i,"class"] = taxon_info.loc[i,"taxon"].split(";")[-3]


#carveme_rank_CF = carveme_rank_CF.sort_values("taxon")
carveme_rank_CF.to_csv(join(model_dir,"rank_carvememodels_CAFBA.csv"))
#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.1'

""" Computing changes in flux distribution in entire intracellular metabolic network
This corresponds to 'theta' in the paper."""

import sys
from pathlib import Path
import argparse
import os
from os.path import join
import glob
from warnings import warn
from multiprocessing import Pool
import numpy as np
import pandas as pd
import cobra

# load CAFBA_mod package
import routine_functions.CAFBA_mod as CF


def cossim_fluxdata(fluxdata):
    '''compute cosine similarity between every set of fluxdata in a dataframe'''
    angle = np.array([])
    for i in range(0,fluxdata.shape[0]):
        flux_i = fluxdata.iloc[i,:].to_numpy()
        for j in range(i+1,fluxdata.shape[0]):
            flux_j = fluxdata.iloc[j,:].to_numpy()
            #j_dist = np.linalg.norm(flux_i-flux_j)
            dotp = flux_i.transpose() @ flux_j
            mod_fluxi = np.power(np.square(flux_i).sum(axis=0), 0.5)
            mod_fluxj = np.power(np.square(flux_j).sum(axis=0), 0.5)
            mods = mod_fluxi*mod_fluxj
            angle = np.insert(angle,len(angle),np.degrees(np.arccos(dotp/mods))) 
    return angle



def cossim_totalflux_multi(target_rxn, NormalizeByGrowth=True):
    
    save_csv = join(savedir,target_rxn + "_fluxchange.csv")
    
    if not os.path.exists(save_csv):
        with open(save_csv, mode = "w") as f:
            f.write(",".join(["exp_id"] + csource_list) + "\n")
    else:
        return
            
    for filepath in glob.glob(join(fluxdir,target_rxn,"*/*" + "_flux.csv")):
        
        exp_id = filepath.split("/")[-1].split("_")[0]
        
        # Flux distribtuion before 
        flux_original = pd.read_csv(glob.glob(join(originaldir,exp_id + "_flux.csv"))[0],
                                    index_col = 0, header = 0)

        # Flux distribution after mutation
        flux_changed = pd.read_csv(filepath,index_col = 0, header = 0)

        # Normalization of the flux distributions by an objective function
        if NormalizeByGrowth:
           flux_changed = flux_changed/flux_changed.loc[obj]
           flux_original = flux_original/flux_original.loc[obj]
        
        flux_change = []
        for s in csource_list:
            if s not in flux_changed.columns.values:
                flux_change.append(np.nan)
                continue
            flux_target_sugar = pd.concat([flux_original[s],flux_changed[s]],axis=1)
            flux_target_sugar = flux_target_sugar.loc[list(set(flux_original.index)&set(flux_changed.index)),:]
            angle = cossim_fluxdata(flux_target_sugar.T)
            flux_change.append(angle[0])
        
        with open(save_csv, mode="a") as f:
            f.write(",".join([exp_id] + [str(x) for x in flux_change]) + "\n")



#################################
#### Main body of the script ####
#################################

if __name__ == "__main__":

    homedir = Path(__file__).parent.absolute()
    
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="directory where the result of FluxSensitivity.py exist")
    p.add_argument("-s","--substrates", type=str, default="None",
                   help="(optional) growth substrates (should be given as a list of 'EX_' reactions by txt file)")
    p.add_argument("-b","--objective",type=str,default="normal", 
                   help="(optional) an  objevtive reaction when computing pfba. \
                   'normal':the iJO1366 biomass funciton, 'no':no normalization,'[reaction id]': The specified reaction is used for the normalization )")
    p.add_argument("-t","--threads", type=int, default=1, help="the number of threads")
    p.add_argument("-o","--optimizer", type=str, default="None", help="(optional) set an optimizer for FBA (e.g., gurobi) Default is cplex.")

    args = p.parse_args()

    # set solver
    if args.optimizer == "None":
        cobra.Configuration().solver = "cplex"
    else:
        cobra.Configuration().solver = args.optimizer

    if args.inputdir== "None":
        print("[WARN] The results directory of randomwalks must be specified by '--inputdir'.")
        print("[WARN] Using test data instead this time.")
        datadir = join(homedir,"test")
    else:
        datadir = args.inputdir
        if not os.path.exists(datadir): 
            sys.stderr.write("%s not exists...\n"%datadir)
            sys.exit()

    
    # set substraets 
    if args.substrates == "None":
        print("[SUBSTRATES] no substrate file is set by -s option.")
        print("[SUBSTRATES] preset 7 sugars will be used.")
        # Carbon sources used in (Takano et al., 2023)
        # I only picked up 7 sugars whose hierarchy matches to experimentally observed one 
        csource_list = ['EX_glc__D_e', 'EX_fru_e', 'EX_man_e','EX_fuc__L_e', 
                      'EX_melib_e','EX_gal_e', 'EX_rib__D_e']
    else:
        substratefile = args.substrates
        with open(substratefile,mode="r") as f:
            csource_list = [x.strip("\n") for x in f.readlines()]

    # load objective function
    if args.objective == "normal":
        # Load original model
        model_e_cafba = CF.CAFBA_mod_Model(
            cobra.io.read_sbml_model(
            join(homedir,"reference models",'iJO1366_cafba.xml')
            )
        )
        obj = [x.id for x in model_e_cafba.reactions if x.objective_coefficient == 1][0]
        NormalizeByGrowth=True
    elif args.objective == "no":
        NormalizeByGrowth=False
    else:
        obj = args.objective
        NormalizeByGrowth=True


    fluxdir = join(datadir,"flux_sensitivity_each")
    originaldir = join(datadir,"flux_original")
    savedir = join(datadir,"flux_cossim_total")
    
    target_set = [x.split("/")[-3] for x in glob.glob(join(fluxdir,"*","*","*" + "_flux.csv"))]
    
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    # Calculate total flux change (cosine similarity)
    with Pool(processes=args.threads) as pool:
        pool.map(cossim_totalflux_multi,target_set)


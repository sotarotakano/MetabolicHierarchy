#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.0'

'''The codes for computing jaccard distance in processing metabolic pathways between 
a given set of carbon sournces. For a given dataset of growth-rate trajectories, the
function computes and returns how the dissmilarity in used metabolic pathways changed
by the additions and deletions of metabolic reactions during randomwalks by reproducing
randomwalks for a given model.'''

import routine_functions.CAFBA_mod as CF
from routine_functions.computing_growth_flux import *
import os
import argparse
from os.path import join
from os.path import basename
import pandas as pd
import numpy as np
import cobra
from multiprocessing import Pool
import sys
from pathlib import Path


def compute_flux_dist(arguments):
    '''Computing the distance in used metabolic reactions (jaccard distance) for a given dataset
    This program is compatible to csv format generated after running "randomwalks.py". '''
    csv_file = arguments[0]
    interval = arguments[1]
    Csource_list = arguments[2]
    basemodel = arguments[3]
    donormodel = arguments[4]
    donormodel2 = arguments[5]
    
    base = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(basemodel))
    donor = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(donormodel))
    donor2 = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(donormodel2))
    #Loading added-deleted reaction list in reference model
    result = pd.read_csv(csv_file,header = 0,index_col = 0)
    result[np.isnan(result)] = 0.0
    
    donor_reaction_id = [x.id for x in donor.reactions]  
    
    rxn_set = result.index
    with base as M:
        for i, rxn in enumerate(rxn_set):
            if "add_" in rxn:
                add_r = rxn.replace("add_","")
                Add_rxns = []
                try:
                    Add_rxns.append(donor.reactions.get_by_id(add_r))
                except:
                    Add_rxns.append(donor2.reactions.get_by_id(add_r))
                if add_r + '_rev' in donor_reaction_id:
                    Add_rxns.append(donor.reactions.get_by_id(add_r + '_rev'))
                M.add_reactions(Add_rxns)

            if "del_" in rxn:
                ko_r = rxn.replace("del_","")
                KO_rxns = []
                KO_rxns.append(M.reactions.get_by_id(ko_r))
                if ko_r + '_rev' in donor_reaction_id:
                    KO_rxns.append(M.reactions.get_by_id(ko_r + '_rev'))
                M.remove_reactions(KO_rxns)
            
            if i%interval == 0:
                flux = computing_flux_pfba(Csource_list,M)
                flux_dist = jaccard_fluxdata(flux.T,min_thresh = 1e-06)
                flux_dist = pd.Series(flux_dist,name = str(i))
                
                if "Flux_dist" not in locals():
                    Flux_dist = flux_dist
                else:
                    Flux_dist = pd.concat([Flux_dist,flux_dist],axis=1)
    
    filename = join(parentdir,"flux_JD_pairwise",csv_file.split("/")[-1].split(".")[0] + '.csv')
    Flux_dist.to_csv(filename)
    print("%s finished."%csv_file)


def jaccard_fluxdata(fluxdata,min_thresh=1e-06):
    '''compute jaccard distance between every set of fluxdata in a dataframe,
    in cafba, forward and reverse reactions are considered as single reaction'''
    Dist = np.array([])
    for i in range(0,fluxdata.shape[0]):
        flux_i = fluxdata.iloc[i,:].to_numpy()
        active_flux_i = {x for x in fluxdata.iloc[i,abs(flux_i) > min_thresh].index}
        active_flux_i = {x.replace("_rev","") for x in active_flux_i}
        for j in range(i+1,fluxdata.shape[0]):
            flux_j = fluxdata.iloc[j,:].to_numpy()
            active_flux_j = {x for x in fluxdata.iloc[j,abs(flux_j) > min_thresh].index}
            active_flux_j = {x.replace("_rev","") for x in active_flux_j}
            j_dist = 1-(len(active_flux_i & active_flux_j)/len(active_flux_i | active_flux_j))
            Dist = np.insert(Dist,len(Dist),j_dist) 
    return Dist


# Main body
if __name__ == "__main__":
    homedir = Path(__file__).parent.absolute()
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-m","--model", type=str, default="None", help="an ancestral model for random walks")
    p.add_argument("-u","--universal",type=str,default="None", 
                   help="an universal model (This will be a pool of reactions for adding)")
    p.add_argument("-r","--resultsdir", type=str, default="None", help="directory where result files of random walks exist")
    p.add_argument("-s","--substrates", type=str, default="None",
                   help="growth substrates (should be given as a list of 'EX_' reactions by txt file)")
    p.add_argument("-t","--threads", type=int, default=1, help="the number of threads")
    p.add_argument("-i","--interval", type=int, default=2000, help="the interval for checking jaccard distance during randomwalks")
    p.add_argument("-o","--optimizer", type=str, default="None", help="(optional) set an optimizer for FBA (e.g., gurobi) Default is cplex.")

    args = p.parse_args()
    # set solver
    if args.optimizer == "None":
        cobra.Configuration().solver = "cplex"
    else:
        cobra.Configuration().solver = args.optimizer

    if args.model== "None":
        print("[MODEL] no model is set by -m option.")
        print("[MODEL] iJO1366 CAFBA is used as an ancestral model.")
        model_e_cafba = join(homedir,'reference models','iJO1366_cafba.xml')
        model_e = join(homedir,'reference models', 'iJO1366.xml')
    else:
        modelfile = args.model
        cafbafile = join(Path(modelfile).parent.absolute(),basename(modelfile).replace(".xml","_cafba.xml"))
        if modelfile.split(".")[-1] == "xml":
            model_e = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(modelfile))
            if not os.path.exists(cafbafile):
                CF.main(modelfile)
            model_e_cafba = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(cafbafile))
        else:
            print("[ERROR] an model file should be given as .xml file format")
            

    if args.substrates == "None":
        print("[SUBSTRATES] no substrate file is set by -s option.")
        print("[SUBSTRATES] preset 7 sugars will be used.")
        # Carbon sources used in (Takano et al., 2023)
        # I only picked up 7 sugars whose hierarchy matches to experimentally observed one 
        Csource_list = ['EX_glc__D_e', 'EX_fru_e', 'EX_man_e','EX_fuc__L_e', 
                      'EX_melib_e','EX_gal_e', 'EX_rib__D_e']
    
    else:
        substratefile = args.substrates
        with open(substratefile,mode="r") as f:
            Csource_list = [x.strip("\n") for x in f.readlines()]

    Npairs = int(len(Csource_list)*(len(Csource_list)-1)/2)
    print("[INPUTS] %i substrates, and %i pairs will be compared."%(len(Csource_list),Npairs))

    if args.universal== "None":
        print("[UNIVERSAL MODEL] no model is set by -m option.")
        print("[UNIVERSAL MODEL] CUiJO1366 is used as a reaction pool.")
        # Universal model (from Bajic et al. PNAS)
        # reversibility and upper and lower bounds are set using 'UiJO1366_default_bounds.csv'
        model_u_e = join(homedir,'reference models', 'CUiJO1366_modified.xml')

        # I made CAFBA model from the above universal model
        model_u_e_cafba = join(homedir,'reference models', 'CUiJO1366_modified_cafba.xml')
    
    else:
        universalfile = args.universal
        ucafbafile = join(Path(universalfile).parent.absolute(),
                          basename(universalfile).replace(".xml","_cafba.xml"))
        if universalfile.split(".")[-1] == "xml":
            model_u_e = universalfile
            if not os.path.exists(ucafbafile):
                CF.main(universalfile)
            model_u_e_cafba = ucafbafile
        else:
            print("[ERROR] an universal file should be given as .xml file format")

    # Setting results directory
    if args.resultsdir== "None":
        print("[WARN] The results directory of randomwalks must be specified by '--resultsdir'.")
        print("[WARN] Using test data instead this time.")
        datadir = join(homedir,"test","results")
    else:
        datadir = args.resultsdir
        if not os.path.exists(datadir): 
            sys.stderr.write("%s not exists...\n"%datadir)
            sys.exit()
    filepath = [join(datadir,x) for x in os.listdir(datadir) if '.csv' in x]
    print("[INPUTS] %i files are found."%len(filepath))

    parentdir = Path(datadir).parent.absolute()

    if not os.path.exists(join(parentdir,"flux_JD_pairwise")):
        os.mkdir(join(parentdir,"flux_JD_pairwise"))
    
    print("[RUNNING] Computing jaccard distance in the processing pathways.")
    arguments = [(x,args.interval,Csource_list,model_e_cafba,model_u_e_cafba,model_u_e) for x in filepath]
    with Pool(processes=args.threads) as pool:
        pool.map(compute_flux_dist, arguments)

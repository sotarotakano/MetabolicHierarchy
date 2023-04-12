#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.1'

'''The python script for flux sensitivity analysis in Takano et al., 2023.
Briefly, we computed the change in flux distribution after adding or deleting
a reaction of interest to a given set of metabolic models.'''

import sys
import argparse
import os
from os.path import join
import glob
from multiprocessing import Pool
import numpy as np
import pandas as pd
import cobra
from cobra.flux_analysis import pfba 
from pathlib import Path
import routine_functions.CAFBA_mod as CF


def computing_growth(EX_C_reactions, model, name = '', lb = -10.0, aerobic = True, normalize_by_Catoms = False, ErrorMessage = False):
    '''compute a growth rate on a given set of carbon sources (EX_C_reactions)'''
    '''return the pandas Series of growth data of a model'''
    '''For optimization, slim_optimize is used.'''
    Growth_matrix = pd.Series(0,index = EX_C_reactions)
    model.set_minimal_media(aerobic = aerobic)
    for i in EX_C_reactions:
        if normalize_by_Catoms:
            try:
                EX_r = model.reactions.get_by_id(i)
                met_formula = next(iter(EX_r.metabolites.keys())).formula
                C_atoms = get_Catoms_number(met_formula)
                if C_atoms != 0:
                    lb_n = lb/C_atoms
                else:
                    lb_n =  lb
            except Exception as exc:
                if ErrorMessage:
                    sys.stderr.write(str(exc) + '\n')
                    lb_n = lb 
        else:
            lb_n = lb
        with model as M:
            try:
                M.reactions.get_by_id(i).lower_bound = lb_n
                growth_rate = M.slim_optimize()
            except Exception as exc:
                if ErrorMessage:
                    sys.stderr.write("Error in slim_optimize ...\n")
                    sys.stderr.write(str(exc) + '\n')  
                growth_rate = 0.0
            if np.isnan(growth_rate):
                growth_rate = 0.0
        Growth_matrix.loc[i] = growth_rate
    Growth_matrix.name = name
    return Growth_matrix


def computing_flux_pfba(EX_C_reactions, model, lb = -20.0, ErrorMessage = True, normalize_by_Catoms = False, internal_reactions = True):
    '''compute flux on a given set of carbon sources (EX_C_reactions)'''
    '''return a pandas dataframe of fluxdata in a given model and a give set of Csources by pfba'''
    model.set_minimal_media()
    
    for i in EX_C_reactions:
        if normalize_by_Catoms:
            try:
                EX_r = model.reactions.get_by_id(i)
                met_formula = next(iter(EX_r.metabolites.keys())).formula
                C_atoms = get_Catoms_number(met_formula)
                if C_atoms != 0:
                    lb_n = lb/C_atoms
                else:
                    lb_n =  lb
            except Exception as exc:
                if ErrorMessage:
                    sys.stderr.write(str(exc) + '\n')
                    lb_n = lb 
        else:
            lb_n = lb
        with model as M:
            if i in M.reactions:
                M.reactions.get_by_id(i).lower_bound = lb_n
                try:
                    fluxdata = pfba(M).fluxes
                    fluxdata.name = i
                except Exception as exc:
                    sys.stderr.write("Error in pfba ...\n")
                    sys.stderr.write(str(exc) + '\n')  
                    continue

        if 'flux_C' not in locals():
            flux_C = fluxdata
        else:
            flux_C = pd.concat([flux_C,fluxdata],axis = 1)           
    
    if 'flux_C' not in locals():
        sys.stderr.write("No flux data is generated from a given set of EX_reactions ...\n")
        return None
    else:
        flux_C = pd.DataFrame(flux_C)
        if internal_reactions:
            rxn_incell = {x.id for x in model.reactions} - model.demand_reactions()
            flux_C  = flux_C.loc[rxn_incell,:] 
        return flux_C


def get_Catoms_number(met_formula):
    Cpos = met_formula.find('C')
    # If there is no C atoms in this molecule, this normalization is skipped,
    # and None is substituted
    if Cpos >= 0:
        num_pos = Cpos
        for k in range(Cpos+1,len(met_formula)):
            if met_formula[k].isnumeric():
                num_pos = k
            else:
                break
        if num_pos == Cpos:
            C_atoms = 1
        else:
            C_atoms = int(met_formula[Cpos+1:num_pos+1]) 
    else:
        C_atoms = 0
    return C_atoms


def flux_original_model(args):
    '''Computing flux distribution in a given set of carbon sources using pFBA.'''
    '''We set lower bounds of EX reactions to -120mM/Catoms'''
    filepath = args[0]
    Sugar_list = args[1]
    
    save_dir = join(datadir,"flux_original")
    exp_id = filepath.split("/")[-1].split("_")[0]

    if os.path.exists(join(save_dir,exp_id + "_flux.csv")):
        return   
    
    # read one of the evolved models
    cafba_mod = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(filepath))
    
    # First computing flux distribution 
    flux = computing_flux_pfba(Sugar_list,cafba_mod,lb = -120,normalize_by_Catoms = True,
                               internal_reactions=False)
    flux.to_csv(join(save_dir,exp_id + "_flux.csv"))


def flux_sensitivity_each(args):
    filepath = args[0]
    target_rxn_set = args[1]
    Sugar_list = args[2]

    save_dir = join(datadir,"flux_sensitivity_each")
        
    exp_id = filepath.split("/")[-1].split("_")[0]
    
    #original_fluxdir =  join(datadir,"flux_original")
    
    # read one of the evolved models
    cafba_mod = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(filepath))
  
    # Add or KO for user defined reactions each. 
    for target_rxn in target_rxn_set:
        error = False
        
        if len(glob.glob(join(save_dir,target_rxn,"*/"+ exp_id + "_flux.csv"))) > 0:
            continue
        
        with cafba_mod as c:
            # If the target reaction exist in the model, the reaction is removed
            if target_rxn in {x.id for x in cafba_mod.reactions}:
                KO_rxns = [cafba_mod.reactions.get_by_id(target_rxn)]
                if target_rxn + '_rev' in {x.id for x in c.reactions}:
                    KO_rxns.append(c.reactions.get_by_id(target_rxn + '_rev'))
                c.remove_reactions(KO_rxns)
                growth = computing_growth(Sugar_list, c, name = Evolved_models[0].split("/")[-1].split(".")[0], 
                                          lb = -120.0, normalize_by_Catoms = True, ErrorMessage = False)
                
                # If this reactions is essential for metabolizing one of a given set of sugars, skip this reaction
                if max(growth) < 1e-06:
                    error = True
                    status = "ESNTL"
                else:
                    after = computing_flux_pfba([x for x in Sugar_list if growth[x] > 1e-06],c,lb = -120,
                                                normalize_by_Catoms = True,internal_reactions=False)
                    status = "KO"
            else:
                # If the target reaction does not exist in the model, the reaction is added
                Add_rxns = [model_u_e_cafba.reactions.get_by_id(target_rxn)]
                if target_rxn + '_rev' in {x.id for x in model_u_e_cafba.reactions}:
                    Add_rxns.append(model_u_e_cafba.reactions.get_by_id(target_rxn + '_rev'))
                c.add_reactions(Add_rxns)
                growth = computing_growth(Sugar_list, c, name = Evolved_models[0].split("/")[-1].split(".")[0], 
                                          lb = -120.0, normalize_by_Catoms = True, ErrorMessage = False)
                # If this reactions is essential for metabolizing one of a given set of sugars, skip this reaction
                if max(growth) < 1e-06:
                    error = True
                    status = "ESNTL"
                else:
                    after = computing_flux_pfba([x for x in Sugar_list if growth[x] > 1e-06],c,lb = -120,
                                                normalize_by_Catoms = True,internal_reactions=False)
                    status = "ADD"
        
        if not error:
            savecsv = join(save_dir, target_rxn, status, exp_id + "_flux.csv")
            if not os.path.exists(join(save_dir,target_rxn)):
                os.mkdir(join(save_dir,target_rxn))

            if not os.path.exists(join(save_dir,target_rxn,status)):
                os.mkdir(join(save_dir,target_rxn,status))

            after.to_csv(savecsv)
    print("Finished %s"%filepath)


def check_flux(args):
    target_rxn = args[0]
    path = args[1]
    NormalizeByGrowth = args[2]
    
    flux = pd.read_csv(path,index_col = 0)
    if NormalizeByGrowth:
        flux = flux/flux.loc[obj]
    
    if target_rxn not in flux.index:
        return [path.split("/")[-1].split("_")[0]]+[np.nan]*len(flux.columns)
    
    if target_rxn + "_rev" in flux.index:
        flux.loc[target_rxn.replace("_rev",""),:] += flux.loc[target_rxn + "_rev",:] 
        flux = flux.drop(target_rxn + "_rev")
    
    target_flux = flux.loc[target_rxn,:]
    
    return [path.split("/")[-1].split("_")[0]]+[x for x in target_flux]



#################################
#### Main body of the script ####
#################################
if __name__ == "__main__":

    homedir = Path(__file__).parent.absolute()
    
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="directory where the result of ModifierScreen.py exist")
    p.add_argument("-s","--substrates", type=str, default="None",
                   help="(optional) growth substrates (should be given as a list of 'EX_' reactions by txt file)")
    p.add_argument("-r","--reactions", type=str, default="None",
                   help="(optional) lists of target reactions (should be given as a list of reaction IDs by txt file). Those reactions are used as a reference (non-modifiers)")
    p.add_argument("-m","--model_dir",type=str,default="None", 
                   help="(optional) directory where the target metabolic models for the analysis exists.")
    p.add_argument("-u","--universal",type=str,default="None", 
                   help="(optional) an universal model (This will be a pool of reactions for adding)")
    p.add_argument("-b","--objective",type=str,default="normal", 
                   help="(optional) an  objevtive reaction when computing pfba. \
                   'normal':the iJO1366 biomass funciton, 'no':no normalization,'[reaction id]': The specified reaction is used for the normalization )")
    p.add_argument("-t","--threads", type=int, default=1, help="the number of threads")
    p.add_argument("-o","--optimizer", type=str, default="None", help="(optional) set an optimizer for FBA (e.g., gurobi) Default is cplex.")
    p.add_argument("-N","--Nmodels", type=int, default=1000, help="The number of models used for the flux sensitivity analysis (default = 1000).")

    args = p.parse_args()

    # set solver
    if args.optimizer == "None":
        cobra.Configuration().solver = "cplex"
    else:
        cobra.Configuration().solver = args.optimizer

    # set target directory
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

    
    # load reference reactions
    if args.reactions == "None":
        print("[REACTIONS] no reaction file is set by -r option.")
        print("[REACTIONS] preset reaction list will be used.")
        #NM_rxnlistfile = join(homedir,"auxiliary","random_NM_list.txt") # full
        NM_rxnlistfile = join(homedir,"auxiliary","random_NM_list_partial.txt") # partial

        with open(NM_rxnlistfile,mode='r') as f:
            non_significant = {x.strip("\n") for x in f.readlines()}
    else: 
        ref_reactionfile = args.reactions
        with open(ref_reactionfile,mode="r") as f:
            non_siginificant = {x.strip("\n") for x in f.readlines()}
    

    # load directory info of target metabolic models
    # Here. we used randomly picked out models generated by randomwalks as a default, but you can specify any directory here.
    if args.model_dir == "None":
        print("[MODELS] no directory is specified by -m option.")
        print("[MODELS] preset metabolic models will be used.")
        model_dir = join(homedir,"test","g_models")
        Evolved_models = sorted([join(model_dir,x) for x in os.listdir(model_dir) if ".xml" in x])
    else: 
        model_dir = args.model_dir
        Evolved_models = sorted([join(model_dir,x) for x in os.listdir(model_dir) if ".xml" in x])
    

    # load universal model (reaction pool)
    if args.universal== "None":
        print("[UNIVERSAL MODEL] no model is set by -m option.")
        print("[UNIVERSAL MODEL] CUiJO1366 is used as a reaction pool.")
        # Universal model (from Bajic et al. PNAS)
        # reversibility and upper and lower bounds are set using 'UiJO1366_default_bounds.csv'
        model_u_e_cafba = CF.CAFBA_mod_Model(
            cobra.io.read_sbml_model(
                join(homedir,'reference models', 'CUiJO1366_modified_cafba.xml')
            )
        )
    else:
        universalfile = args.universal
        if os.path.exists(universalfile):
            if universalfile.split(".")[-1] == "xml":
              model_u_e_cafba = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(universalfile))
            else:
                print("[ERROR] an universal file should be given as .xml file format")
        else:
            print("[ERROR] %s does not exist. can't load the model"%universalfile)
    

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



    ''' Loading capacitor potentiator data'''

    '''We focused on the flux sensitivity against the mutations on evolutionary modifiers,
    we first listed those reactions by referring "CP_individual_flip.csv", which is the
    output of "ModifierScreen.py". For comparing the impacts of those modifiers to those of
    non-modifiers, we randomly picked out 1000 reactions from non-modifiers. For securing the
    reproducibility, we used preliminary selected datasets ("random_NM_list.txt").'''

    CP_individual_fl = pd.read_csv(join(datadir,"Summary_modifiers_individual.csv"), index_col = 0)
    met_list = [x for x in CP_individual_fl.columns.values if 
                "EX_" + x in csource_list and "pval" not in x and x != "Modifier"]
    
    if len(met_list) < 1:
        print("No substrates matched to the result of 'Summary_modifiers_individual.csv'.")
        sys.exit()
    
    pvalue_mets = ["pval_" + m for m in met_list]

    significant_fl = CP_individual_fl.loc[[x for x in CP_individual_fl.index 
                                        if min(CP_individual_fl.loc[x,pvalue_mets]) < 0.01 
                                        and min(abs(CP_individual_fl.loc[x,met_list])) > 1e-06],:]

    significant = {x for x in significant_fl.index}

    # non_significant = {x for x in CP_individual_fl.index if x not in significant}  # optional

    target_set = significant | non_significant

    print("Target reactions are %i reactions"%len(target_set))


    ### Computing flux distribution ###
    print("[RUNNING] Computing flux sensitivity...")

    if not os.path.exists(join(datadir,"flux_original")):
        os.mkdir(join(datadir,"flux_original"))
    
    ''' Calculate flux before mutations'''
    arguments = ([e,csource_list] for e in Evolved_models[:min(len(Evolved_models),args.Nmodels)])
    with Pool(processes=args.threads) as pool:
        pool.map(flux_original_model, arguments) 


    # Check flux magnitude of reaction
    print("[RUNNING] Checking flux magnitude...")
    savecsv = join(datadir,"flux_magnitude.csv")
    columns = ["reaction"] + csource_list
    with open(savecsv, mode="w") as csv:
        csv.write(",".join(columns)+"\n")

    filelist = glob.glob(join(datadir,"flux_original","*" + "_flux.csv"))
    print("There are %i files of flux distributions after random walks"%len(filelist))
    for rxn in target_set:
        with Pool(processes=args.threads) as pool:
            flux = pool.map(check_flux, 
                            [[rxn,f, NormalizeByGrowth] for f in filelist])
        flux = pd.DataFrame(flux, columns = ["exp_id"] + csource_list).set_index("exp_id").dropna()
        flux_mean = flux.mean()
        with open(savecsv, mode="a") as csv:
            csv.write(",".join([rxn] + [str(x) for x in flux_mean])+"\n")

    
    ## Analyzing flux distribution after additon or deletion ##
    if not os.path.exists(join(datadir,"flux_sensitivity_each")):
        os.mkdir(join(datadir,"flux_sensitivity_each")) 

    ''' Calculate flux after mutations '''
    arguments = ([e,target_set, csource_list] for e in Evolved_models[:min(len(Evolved_models),args.Nmodels)])
    with Pool(processes=args.threads) as pool:
        pool.map(flux_sensitivity_each, arguments) 

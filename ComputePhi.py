#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.1'

""" Computing changes in flux distribution on a specific metabolites. 
Briefly, the program calculated cosine similarity in a influx-efflx matrix 
on an individual metabolite between before and after a mutation on a target reaction. 
This metric corresponds to 'phi' in Takano et al., 2023 biorxiv. """

import os
import sys
import argparse
from pathlib import Path
from os.path import join
from multiprocessing import Pool
import numpy as np
import pandas as pd
import cobra
import glob
import routine_functions.CAFBA_mod as CF


def cossim_met(model,metid,flux_wt,flux_mut,direction=False):
    
    '''Computing cosine similarity in the flux on the target metablite 
    between the given flux distributions.'''

    if not set(flux_wt.columns.values) == set(flux_mut.columns.values):
        print("Two flux matrices have different columns.")
        return
    
    # Picking out reactions which can be involved in the target metabolite
    met_rxn = {x.id for x in list(model.metabolites.get_by_id(metid).reactions)}
    met_rxn = list(met_rxn & {x for x in flux_wt.index} & {x for x in flux_mut.index})
    met_flux_wt = flux_wt.loc[met_rxn,:]
    met_flux_mut = flux_mut.loc[met_rxn,:]
    
    
    # (Optional) Making the influx (i.e., production of a target metabolite) is positive
    # and efflux (i.e., comsumption of a target metabolite) is negative
    if direction:
        met = model.metabolites.get_by_id(metid)
        for r in met_rxn:
            mets = model.reactions.get_by_id(r).metabolites
            try:
                direction = mets[met]
                if direction > 0:
                    pass
                else:
                    met_flux_wt.loc[r,:] *= -1
                    met_flux_mut.loc[r,:] *= -1
            except Exception as exc:
                print(exc)
                continue

    # Computing cosine similarity 
    angle = []
    for j in flux_wt.columns.values:
        muts = np.array(met_flux_mut.loc[:,j])
        wt = np.array(met_flux_wt.loc[:,j])
        dotp = muts.transpose() @ wt
        mod_muts = np.power(np.square(muts).sum(axis=0), 0.5)
        mod_wt = np.power(np.square(wt).sum(axis=0), 0.5)
        mods = mod_muts*mod_wt

        if mods < 1e-12:
            angle.append(np.nan)
        elif dotp/mods > 1:
            angle.append(np.nan)
        else:
            angle.append(np.degrees(np.arccos(dotp/mods)))
    
    angle = np.array(angle)
    angle[np.isnan(angle)] = 0.0  
    
    return angle


def calculate_phi(args):

    ''' Calculating the impact on flux distribution on individual metabolites
    when adding or deleting a given target reaction. The changes in flux distribution 
    between before and after the mutations are computed using all existing files.'''

    target_rxn = args[0]
    NormalizeByGrowth = args[1]

    if not os.path.exists(join(savedir,target_rxn)):
        os.mkdir(join(savedir,target_rxn))  
    
    # The files of flux distributions after addition or deletion of the target reaction
    flux_mut_files = glob.glob(join(fluxdir,target_rxn,"*/*" + "_flux.csv"))

    for filepath in flux_mut_files:
        exp_id = filepath.split("/")[-1].split("_")[0]
        
        status = filepath.split("/")[-2]
        # Original flux distribution
        flux_original = pd.read_csv(glob.glob(join(originaldir,exp_id + "_flux.csv"))[0],
                                    index_col = 0, header = 0)
        flux_changed = pd.read_csv(filepath,index_col = 0, header = 0)
        
        if not set(flux_original.columns.values) == set(flux_changed.columns.values):
            #print("Two flux matrices have different columns. expid = %s; reaction = %s"%(exp_id,target_rxn))
            continue

        # Normalized by growth rate
        if NormalizeByGrowth:
            flux_changed = flux_changed/flux_changed.loc[obj]
            flux_original = flux_original/flux_original.loc[obj]

        # To reduce the calculation cost, only the changed flux is focused.
        delta_flux = flux_changed - flux_original
        delta_flux = delta_flux.loc[abs(delta_flux).mean(axis=1) > 1e-06,:]
        significant = delta_flux.index
        
        # Then we only focus on the metabolites involved in the reactions whose flux were changed
        significant_mets = set()
        
        for s in significant:
            significant_mets = significant_mets | {x.id for x in model_u_e_cafba.reactions.get_by_id(s).metabolites}

    
        # Compute the flux-impact for the active metabolites
        for metid in significant_mets:

            save_csv = join(savedir,target_rxn,metid + "_fluxchange.csv")
            if not os.path.exists(save_csv):
                with open(save_csv, mode = "w") as f:
                    f.write(",".join(["exp_id","status"] + csource_list) + "\n")

            with open(save_csv,mode="r") as f:
                saveddata = f.readlines()
                if len([x for x in saveddata if exp_id in x]) > 0:
                    continue

            flux_change = cossim_met(model_u_e_cafba, metid, flux_original, flux_changed, direction=False)
         
            with open(save_csv, mode="a") as f:
                f.write(",".join([exp_id,status] + [str(x) for x in flux_change]) + "\n")

        
def comparison_phi(args):

    met = args[0]
    status = args[1]
    datadir = args[2]
    savedir = args[3]
    
    filelist = glob.glob(join(datadir,"*" + "metsummary.csv"))

    savecsv = join(savedir,met + "-"+ status +"mean.csv")
    if not os.path.exists(savecsv):
        with open(savecsv, mode = "w") as f:
            f.write(",".join(["reaction","sugar","value"]) + "\n")
    else:
        #print("%s already exists"%savecsv)
        return
            
    for f in filelist:
        rxn = f.split("/")[-1].replace("_metsummary.csv","")
        fc_dt = pd.read_csv(f,index_col = 0, header = 0)
        target = fc_dt.loc[fc_dt["metabolite"] == met,:]
        target = target.loc[target["C-P"] == status,:]
        for t in target.index:
            with open(savecsv, mode="a") as f:
                factor = summary_totalFX.loc[(summary_totalFX["sugar"]==target.loc[t,"sugar"])
                                             *(summary_totalFX["reaction"]==target.loc[t,"reaction"]),"fx_cossim"]
                nom_fxsens = target.loc[t,"fx_mean"]/factor.iloc[0]
                f.write(",".join([rxn,target.loc[t,"sugar"],str(nom_fxsens)])+"\n")



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
    p.add_argument("-u","--universal",type=str,default="None", 
                   help="an universal model. This is necessary to access metabolite info for each reaction. Default is CUiJO1366 (same as our paper)")
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
    
    # Load universe model
    if args.universal== "None":
        print("[UNIVERSAL MODEL] no model is set by -m option.")
        print("[UNIVERSAL MODEL] CUiJO1366 is used as a reaction pool.")
        # Universal model (from Bajic et al. PNAS)
        # reversibility and upper and lower bounds are set using 'UiJO1366_default_bounds.csv'
        model_u_e = CF.CAFBA_mod_Model(
            cobra.io.read_sbml_model(
                join(homedir,'reference models', 'CUiJO1366_modified.xml')
            )
        )

        # I made CAFBA model from the above universal model
        model_u_e_cafba = CF.CAFBA_mod_Model(
            cobra.io.read_sbml_model(
                join(homedir,'reference models', 'CUiJO1366_modified_cafba.xml')
            )
        )
        model_u_e_cafba.add_reactions([model_u_e.reactions.get_by_id("G3PD1")])
        model_u_e_cafba.reactions.get_by_id("G3PD1")

    else:
        universalfile = args.universal
        if universalfile.split(".")[-1] == "xml":
            if "cafba" not in universalfile:
                print("[WARN] This file does not include 'cafba' in the file name.")
            model_u_e_cafba = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(universalfile))
        else:
            print("[ERROR] an universal file should be given as .xml file format")


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


    originaldir = join(datadir,"flux_original")
    fluxdir = join(datadir,"flux_sensitivity_each")

    # target reactions: Here, the program will targets the reactions whose flux distribution were computed 
    # and saved in the directory, "flux_sensitivity_each".
    target_set = [x.split("/")[-3] for x in glob.glob(join(fluxdir,"*","*","*" + "_flux.csv"))]

    print("Target reactions are %i reactions"%len(target_set))

    savedir = join(datadir,"flux_cossim_mets")
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    with Pool(processes=args.threads) as pool:
        pool.map(calculate_phi, [(x,NormalizeByGrowth) for x in target_set])
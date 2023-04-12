
from os.path import join
import os
import gc
import glob
from scipy import stats
import statsmodels.stats.multitest
from multiprocessing import Pool
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import sys


def phi_by_reaction(args):

    '''Merging the data of flux sensitivity score (phi) by reactions.
    This fucntion computes the average phi score of a target reaction in
    all examined genetic backgrounds across the metabolites.'''

    target_rxn = args[0]
    target_mets = args[1]
    CP_individual = args[2]
    
    if os.path.exists(join(savedir,target_rxn + "_metsummary.csv")):
        print("%s already finished"%target_rxn)
        return

    fc_dt = []
    for m in target_mets:
        for j in cmets_list:
            s = "EX_" + j
            
            if target_rxn in CP_individual.index:
                pvalue = CP_individual.loc[target_rxn,"pval_"+ j]
                value = CP_individual.loc[target_rxn,j] 

                if pvalue < 0.01:
                    if value > 0:
                        status = "C"
                    else:
                        status = "P"
                else:
                    status = "non_CP"
            else:
                value = np.nan
                status = "Unknown"

            
            if os.path.exists(join(fluxmetdir,target_rxn, m + "_fluxchange.csv")):
                sensitivity_df = pd.read_csv(join(fluxmetdir,target_rxn,m + "_fluxchange.csv"), 
                                                index_col = 0, header = 0) 
            else:
                #print("not exist %s"%join(savedir,target_rxn, m + "_fluxchange.csv"))
                continue
            sensitivity_df = sensitivity_df.sort_index()


            fx_mean = sensitivity_df.loc[:,s].mean()  # the average phi in the all genetic backgrounds
            newdata = [target_rxn,m,j,value,status,fx_mean,
                        len(sensitivity_df.loc[~np.isnan(sensitivity_df[s])])]

            fc_dt.append(newdata)

    fc_dt = pd.DataFrame(fc_dt,columns=["reaction","metabolite","sugar","delta_theta","C-P",
                                        "fx_mean","appearance"])
    fc_dt = fc_dt.loc[~np.isnan([x for x in fc_dt["fx_mean"]]),:]
    
    if len(fc_dt) > 0:
        fc_dt.to_csv(join(savedir,target_rxn + "_metsummary.csv"))
    else:
        pass


def comparison_phi(args):

    '''Summarizing phi score by metabolites. One of our objective is screening
    significantly sensitive metabolites specifically by the mutations on modifiers (C-P reactions)
    compared to non modifiers. Here, this code returns the phi score of a target metabolite against
    mutations on either of "C"(Capacitor), "P"(Potentiator), or "non_CP" (non modifiers).
    Users can specify that by changing 'status' variable. In default, the codes below compute
    all three cases, enabling the compatison among all three cases. '''

    met = args[0]
    status = args[1]
    savedir = args[2]
    
    filelist = glob.glob(join(savedir,"*" + "metsummary.csv"))

    if len(filelist) < 1:
        print("[WARN] no summaryfile exists...")
        return

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
            with open(savecsv, mode="a") as csv:
                csv.write(",".join([rxn,target.loc[t,"sugar"],str(target.loc[t,"fx_mean"])])+"\n")


# Main body 
if __name__ == "__main__":

    homedir = Path(__file__).parent.absolute()
    
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="directory where the result of ModifierScreen.py exist")
    p.add_argument("-s","--substrates", type=str, default="None",
                   help="(optional) growth substrates (should be given as a list of 'EX_' reactions by txt file)")
    p.add_argument("-t","--threads", type=int, default=1, help="the number of threads")

    args = p.parse_args()


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


    # target reactions: Here, I selected reactions which have significantly high delta theta value (C-P reactions)
    # and randomly selected other reactions.
    fluxdir = join(datadir,"flux_sensitivity_each")
    target_set = {x.split("/")[-3] for x in glob.glob(join(fluxdir,"*","*","*" + "_flux.csv"))}
    

    savedir= join(datadir,"phi_summary")
    if not os.path.exists(savedir):
        os.mkdir(savedir)

    ''' Loading capacitor potentiator data'''

    '''We focused on the flux sensitivity against the mutations on evolutionary modifiers,
    we first listed those reactions by referring "CP_individual_flip.csv", which is the
    output of "ModifierScreen.py". For comparing the impacts of those modifiers to those of
    non-modifiers, we randomly picked out 1000 reactions from non-modifiers. For securing the
    reproducibility, we used preliminary selected datasets ("random_NM_list.txt").'''

    CP_individual_fl = pd.read_csv(join(datadir,"Summary_modifiers_individual.csv"), index_col = 0)
    cmets_list = [x for x in CP_individual_fl.columns.values if 
                "EX_" + x in csource_list and "pval" not in x and x != "Modifier"]

    if len(cmets_list) < 1:
        print("No substrates matched to the result of 'Summary_modifiers_individual.csv'.")
        sys.exit()

    pvalue_mets = ["pval_" + m for m in cmets_list]

    significant_fl = CP_individual_fl.loc[[x for x in CP_individual_fl.index 
                                        if min(CP_individual_fl.loc[x,pvalue_mets]) < 0.01 
                                        and min(abs(CP_individual_fl.loc[x,cmets_list])) > 1e-06],:]

    significant = {x for x in significant_fl.index}


    ## Summarizing cosine similarity results

    # All metabolites
    fluxmetdir = join(datadir,"flux_cossim_mets")
    mets_list = {x.split("/")[-1].replace("_fluxchange.csv","") for x 
                in glob.glob(join(fluxmetdir,"*", "*"+ "_fluxchange.csv"))}


    print("Summarizing flux impact on metabolites in individual reaction")

    # Create summaries by reactions
    arguments = [(x,mets_list,CP_individual_fl) for x in target_set]
    with Pool(processes=args.threads) as pool:
        pool.map(phi_by_reaction,arguments)


    # Create summaries by metabolites
    with Pool(processes=args.threads) as pool:
        pool.map(comparison_phi,[(x,"C",savedir) for x in mets_list])

    with Pool(processes=args.threads) as pool:
        pool.map(comparison_phi,[(x,"P",savedir) for x in mets_list])

    with Pool(processes=args.threads) as pool:
        pool.map(comparison_phi,[(x,"non_CP",savedir) for x in mets_list])


    '''Screening sensitive metabolites'''
    '''The following codes statistically screen the significantly sensitive mets
    against the mutation on modifiers.'''
    
    met_pvalue = pd.DataFrame(np.nan,index=mets_list,columns=["C","P","CP","nonCP","CP vs N","C vs P"])
    target_mets = {x.split("/")[-1].split("-")[0] for x 
                   in glob.glob(join(savedir,"*" + "-non_CPmean.csv"))}
    
    target_mets = {x for x in target_mets if os.path.exists(join(savedir,x + "-Cmean.csv"))
                   and os.path.exists(join(savedir,x + "-Pmean.csv")) and "_c" in x}
    
    min_thresh = 1e-06

    for t in target_mets:
        sensitivity_C = pd.read_csv(join(savedir,t + "-Cmean.csv"), header = 0)
        sensitivity_P = pd.read_csv(join(savedir,t + "-Pmean.csv"), header = 0)
        sensitivity_non = pd.read_csv(join(savedir,t + "-non_CPmean.csv"), header = 0)
            
        met_pvalue.loc[t,"C"] = max(np.mean(sensitivity_C["value"]),min_thresh)
        met_pvalue.loc[t,"P"] = max(np.mean(sensitivity_P["value"]),min_thresh)
        met_pvalue.loc[t,"CP"] = max(np.mean(pd.concat([sensitivity_C["value"],sensitivity_P["value"]])),min_thresh)
        met_pvalue.loc[t,"nonCP"] = max(np.mean(sensitivity_non["value"]),min_thresh)
        
        if len(sensitivity_C) > 1 and len (sensitivity_P) > 1:
            met_pvalue.loc[t,"C vs P"] = stats.ranksums(sensitivity_C["value"],sensitivity_non["value"])[1]
        
        if len(sensitivity_non) > 1: 
            if len(sensitivity_P) > 1 or len(sensitivity_C) > 1:
                met_pvalue.loc[t,"CP vs N"] = stats.ranksums(pd.concat([sensitivity_C["value"],sensitivity_P["value"]]),
                                                            sensitivity_non["value"])[1]

        
    met_pvalue = met_pvalue.dropna(how="any",axis = 0)
    # p-value correction
    cor_pval_CP = statsmodels.stats.multitest.multipletests(met_pvalue["C vs P"],method="hs")[1]
    cor_pval_All = statsmodels.stats.multitest.multipletests(met_pvalue["CP vs N"],method="hs")[1]
    met_pvalue["C vs P"] = cor_pval_CP
    met_pvalue["CP vs N"] = cor_pval_All

    # Calculate Difference
    met_pvalue["diff(C-P)"] = met_pvalue["C"] - met_pvalue["P"]
    met_pvalue["diff(CP-N)"] = met_pvalue["CP"] - met_pvalue["nonCP"]

    met_pvalue.to_csv(join(datadir,"Sensitive_mets.csv"))
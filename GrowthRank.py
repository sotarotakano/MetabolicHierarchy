#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.0'

'''The codes for converting growth rate matrix to the growth rate rank 
for a given set of carbon sources during randomwalks. Before running this program,
"randomwalks.py" should be run for getting the set of evolutionary trajectories.'''

import copy
import os
import argparse
from os.path import join
import pandas as pd
import numpy as np
import sys
from multiprocessing import Pool
from warnings import warn
from pathlib import Path
import pickle
import gc


###>>>> Rank computing
def compute_rankdata(arguments):
    filepath = arguments[0]
    ow = arguments[1]
    non_zero = arguments[2]

    '''Return the ranks in the growth rate for a given set of growth rate
    matrix during random walks.'''
    result = pd.read_csv(filepath, header = 0, index_col = 0)
    result[np.isnan(result)] = 0.0
    savefile = join(parentdir,"rankresult","rank_"+filepath.split("/")[-1].split(".")[0]+".csv")
    if os.path.exists(savefile) and not ow:
        return
    if min(result.iloc[-1,:]) > 1e-06:
        rankdata = rank_data_matrix(result,min_threshold = 1e-03,non_zero = non_zero)
        rankdata.to_csv(savefile)


def rank_data_matrix(dataframe,min_threshold = 1e-05,non_zero = False):
    if non_zero: # whether non-growth-supporting sugar is included in a rank
        dataframe = dataframe.loc[:,dataframe.iloc[-1,:] > 0]
    Carbon_rank = copy.copy(dataframe)
    for N in range(0,dataframe.shape[0]):
        growth_data = copy.copy(dataframe.iloc[N,:]).to_numpy()
        growth_data[growth_data < 1e-06] = 1e-06
        rank_data = np.zeros(len(growth_data))
        i, k = 1, 0
        while min(growth_data) < 10000:
            min_growth = min(growth_data)
            rank_data[growth_data == min_growth] = i
            k += len(growth_data[growth_data == min_growth])
            growth_data[growth_data == min_growth] = 10000
            if np.log10(min(growth_data)/min_growth) > min_threshold:
                i = i + k
                k = 0
        Carbon_rank.iloc[N,:] = rank_data
    return Carbon_rank
###<<<<< END rank computing


###>>>>> Computing rank shifts
def rank_shift(filepath):
    '''Count the total number of rank shifts (moving distance) of individual Csources
    in each randomwalk-experiment. If the Csource of interest move from "2" to "4" after a
    certain mutation event, that is regarded as 2 rank-shifts. This function returns such 
    rank shifts in total during randomwalks.'''

    try:
        rankresult = pd.read_csv(join(datadir,filepath),header = 0,index_col =0)
    except Exception as exc:
        sys.stderr.write("Error in load data from csvfile ...\n")
        sys.stderr.write(str(exc) + '\n')
        return
    rankresult = rankresult.T
    EX_list = rankresult.index
        
    data = []
    for rxn in EX_list:
        # rank result data is converted to numpy
        rank_timecourse = rankresult.loc[rxn,:].to_numpy()
        rank_change_timecourse = abs(np.diff(rank_timecourse))
        data.append(np.sum(rank_change_timecourse))
            
    data = pd.Series(data = data,index = EX_list,name = filepath.split("/")[-1].split(".")[0].replace("rank_",""))
    return data


def rank_shift_multi(datadir,savecsv,threads):
    filepath = [join(datadir,x) for x in os.listdir(datadir) if '.csv' in x]
    
    with Pool(processes=threads) as pool:
        data_all = pool.map(rank_shift, filepath)
    
    print("Calculation finished...")
    
    for i in range(len(data_all)):
        if i == 0:
            with open(savecsv, mode="w") as f:
                f.write(","+",".join([x for x in data_all[i].index]) + "\n")
        
        with open(savecsv, mode="a") as f:
            f.write(data_all[i].name + "," + ",".join([str(x) for x in data_all[i]]) + "\n")
        
        if i%100 == 0 and i > 0:
            print("%i trajectories finished"%i)
        
    return

###<<<<< END Computing rank shifts


###>>>>> Computing rank filps (pairwise)
def pairwise_flipcount_multi(filepath):
    '''Compute how many rank flip events occur among all pairs of substrates in a dataframe.
    The result is saved in a  .pkl format'''
    save_dir = join(parentdir,'flipcount_pairwise')
    savefile = join(save_dir,'rankflip_' + filepath.split("/")[-1].split(".")[0].split("_")[-1] + '.pkl')
    
    if os.path.exists(savefile):
        return
    
    try:
        rankresult = pd.read_csv(filepath,header = 0,index_col =0)
        rankresult[np.isnan(rankresult)] = 0.0
          
    except Exception as exc:
        sys.stderr.write("Error in load data from csvfile ...\n")
        sys.stderr.write(str(exc) + '\n')  
        return
    
    rankdata = rankresult
    flip_array = flip_count_array(rankdata)
    
    filehandler = open(savefile, 'wb')
    pickle.dump(np.array(flip_array), filehandler)
    filehandler.close()
    gc.collect()
    return np.array(flip_array)


def flip_count_array(rankdata):
    '''compute how many rank flip events occur among all pairs of substrates in a dataframe
    and return as a numpy array'''
    fliptotal = []
    for i in range(0,rankdata.shape[1]):
        rank_i = rankdata.iloc[:,i].to_numpy()
        for j in range(i+1,rankdata.shape[1]):
            rank_j = rankdata.iloc[:,j].to_numpy()
            flip = flip_count(rank_i,rank_j)
            fliptotal.append(flip)
    fliptotal = np.array(fliptotal)
    return fliptotal


def flip_count(rankdata1,rankdata2):
    '''Counting the number of flip events between the rank of two carbon sources in deletion experiments'''
    if len(rankdata1) == len(rankdata2):
        preference = np.zeros(len(rankdata1))
        preference[np.where(rankdata1 < rankdata2)[0]] = 1
        preference[np.where(rankdata1 > rankdata2)[0]] = -1
        diff = np.diff(preference)
        flip = len(diff[diff != 0.0])
        return flip
    else:
        raise Exception('rankdata1 and rankdata2 should be same length.' +
                        'rankdata1 = %i, rankdata2 = %i ' % (len(rankdata1),len(rankdata2)))
    
###<<<<< END Computing rank flips (pairwise)

# Main body
if __name__ == "__main__":

    homedir = Path(__file__).parent.absolute()
    
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-i","--inputdir", type=str, default="None", help="directory where result files of random walks exist")
    p.add_argument("-t","--threads", type=int, default=1, help="the number of threads")
    p.add_argument("-ow","--over_write", action='store_true')
    p.add_argument("-gs","--growth_substrate", action='store_true',help="focusing on only substrates supporting the growth by the end of random walks")

    args = p.parse_args()

    if args.over_write:
        ow = args.over_write 
    else:
        ow = False
    
    if args.growth_substrate:
        non_zero = args.growth_substrate 
    else:
        non_zero = False

        
    if args.inputdir== "None":
        print("[WARN] The results directory of randomwalks must be specified by '--inputdir'.")
        print("[WARN] Using test data instead this time.")
        datadir = join(homedir,"test","results")
    else:
        datadir = args.inputdir
        if not os.path.exists(datadir): 
            sys.stderr.write("%s not exists...\n"%datadir)
            sys.exit()

    filepath = [join(datadir,x) for x in os.listdir(datadir) if '.csv' in x]
    if len(filepath) < 1:
        warn("no files exist in %s."%datadir)
        sys.exit()
    
    print("[INPUTS] %i files are found."%len(filepath))

    parentdir = Path(datadir).parent.absolute()

    if not os.path.exists(join(parentdir,"rankresult")):
        os.makedirs(join(parentdir,"rankresult"))

    print("[RUNNING] Computing growth-rate ranks for %s files"%len(filepath))

    # Convert growth rate data to growth rank data
    arguments = [(x,ow,non_zero) for x in filepath]
    with Pool(processes=args.threads) as pool:
        pool.map(compute_rankdata, arguments)

    print("[RUNNING] Computing rank shifts...")

    # Computing rank changes
    rank_shift_multi(join(parentdir,"rankresult"),
                     join(parentdir,"rankchange_individual_sugars.csv"),
                     args.threads)

    if not os.path.exists(join(parentdir,"flipcount_pairwise")):
        os.makedirs(join(parentdir,"flipcount_pairwise"))

    print("[RUNNING] Computing rank flips in pairwise manner...")

    # Computing rank filps (pairwise manner)
    with Pool(processes=args.threads) as pool:
        pool.map(pairwise_flipcount_multi, filepath)

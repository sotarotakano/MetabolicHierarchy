#!/usr/bin/env python
__author__ = 'Sotaro Takano and Djordje Bajic'
__version__ = '1.0.2'

'''The codes for generating randomly mutated models (Takano et al., 2023, biorxiv). 
This program add and delete metabolic reactions to the models under the condition 
that they can grow on a given set of carbon sources.
All the codes are compatible with cobrapy and CAFBA model.'''


import routine_functions.CAFBA_mod as CF
import sys
import os
from os.path import join
from os.path import basename
import pandas as pd
import numpy as np
import cobra
from multiprocessing import Pool
import secrets
import datetime
import argparse
from pathlib import Path

# Functions used for computing growth

def computing_growth(EX_C_reactions, model, name = '', lb = -10.0,
                     aerobic = True, normalize_by_Catoms = False, ErrorMessage = False):
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


# Functions used for deletion experiments

def random_add(EX_reactions, model, Add_rxns, min_threshold,
                normalize_by_Catoms=False, lb=-10.0):
    ''' Check whether a random addition event satisfies given constraints
        : growth on a given set of Csources (normally)
        Then, returns the appendability and list of growth rate on those Csources'''

    # From Bajic et al., PNAS. 
    # Here, addition-deletion experiments are performed under the conditon that
    # any of the following reaction pairs don't coexist in the model.
    # For instance, if SHSL2 already exists in the model, 
    # an addition of SHSL2r is not allowed in the experiments.
    excluded_pairs = [('SHSL2','SHSL2r'),('DHORD_NAD','DHORDi'),
    ('ENO','HADPCOADH'),('LEUTA','LLEUDr'),('P5Crx','PRO1y')]

    with model as M:
        try:
            M.add_reactions(Add_rxns)
        except Exception as exc:
            sys.stderr.write("Error in knock out reaction ...\n")
            sys.stderr.write(str(exc) + '\n')
            return False, None, None
        rxn_list = [x.id for x in M.reactions]
        non_appendable = [x for x in excluded_pairs if x[0] in rxn_list and x[1] in rxn_list]
        if len(non_appendable) > 0:
            print("%s and %s coexist... this addition event canceled..." % (non_appendable[0][0],non_appendable[0][1]))
            return False, None, None
        else:
            growth_data = computing_growth(EX_reactions, M, '', normalize_by_Catoms=normalize_by_Catoms,lb=lb)
            Avail_C = len(growth_data[growth_data > 1e-06])
            if Avail_C > min_threshold:
                return True, growth_data, Avail_C
            else:
                return False, growth_data, Avail_C


def random_KO(EX_reactions, model, KO_rxns, min_threshold,
                normalize_by_Catoms=False, lb=-10.0):
    ''' Check whether a random deletion event satisfies given constraints
        : growth on  a given set of Csources (normally)
        Then, returns the deletability and list of growth rate on those Csources'''

    with model as M:
        try:
            M.remove_reactions(KO_rxns)
        except Exception as exc:
            sys.stderr.write("Error in knock out reaction ...\n")
            sys.stderr.write(str(exc) + '\n')
            return False, None, None
        growth_data = computing_growth(EX_reactions, M, '', normalize_by_Catoms=normalize_by_Catoms,lb=lb)
        Avail_C = len(growth_data[growth_data > 1e-06])
        if Avail_C > min_threshold:
            return True, growth_data, Avail_C
        else:
            return False, growth_data, Avail_C


def long_term_addition_deletion(model, donor, EX_reactions, possible_dels, possible_adds, Ntrials,
                                min_threshold,save_folder,filename,lb=-10.0, RunToFinal=True, 
                                normalize_by_Catoms=False):
    '''Main body of the long-term deletion and HGTexperiment algorithm.
       Results of the experiment are written in csv file'''
    
    csvname = save_folder + '/' + filename + '.csv'
    digit = int(np.sqrt(sum([ord(x)**n for n,x in enumerate(filename)])))
    
    #List of deletable reactions in a model
    possible_adds = list(possible_adds)
    possible_dels = list(possible_dels)
    
    # id list of possible additional reactions in a donor(universal) model
    donor_reaction_id = [x.id for x in donor.reactions]    

    trial = 0

    with model as M:
    #Compute growth in a original model   
        growth_data = computing_growth(EX_reactions,M,'anc')
        EX_list = growth_data.index.to_list()
        
        if not os.path.exists(csvname):
            with open(csvname,mode='wt') as nf:
                nf.write(',' + ','.join(EX_list) + '\n')

        print('Addition-Deletion experiment started...')
        
        while trial < Ntrials:
            # id list of reactions in the ancestral model
            model_reaction_id = [x.id for x in M.reactions]

            # Deleting randomly chosen reactions in the model
            done = 0
            search = 0

            while not done:
                np.random.seed(digit+trial+search)
                KO_rxns = []
                # Get one candidate of knock-out reaction (and its reverse reaction) randomly from original ecoli model
                ko_r = possible_dels[np.random.randint(0, len(possible_dels))]
                KO_rxns.append(M.reactions.get_by_id(ko_r))
                if ko_r + '_rev' in model_reaction_id:
                    KO_rxns.append(M.reactions.get_by_id(ko_r + '_rev'))
                done, growth_data, Avail_C = random_KO(EX_reactions, M, KO_rxns,  
                min_threshold, normalize_by_Catoms=normalize_by_Catoms, lb=lb)
                search += 1
           
            # After confirming that knocking out the randomly chosen reaction satisfies constraints,
            # this reaction is truly deleted from the model.
            
            M.remove_reactions(KO_rxns)
            
            # Changing possible reaction list for addition and deletion
            possible_dels.remove(ko_r)
            possible_adds.append(ko_r)
            ko_result = ["del_" + ko_r] + [str(x) for x in growth_data]

            if os.path.exists(csvname):
                with open(csvname,mode='a') as f:
                    f.write(','.join(ko_result) + '\n')


            # Adding randomly chosen reactions to the model
            done = 0
            while not done:
                Add_rxns = []
                # Get one reaction (and its reverse reaction) randomly from universal reaction pool
                add_r = possible_adds[np.random.randint(0,len(possible_adds))]
                Add_rxns.append(donor.reactions.get_by_id(add_r))
                if add_r + '_rev' in donor_reaction_id:
                    Add_rxns.append(donor.reactions.get_by_id(add_r + '_rev'))

                done, growth_data, Avail_C = random_add(EX_reactions, M, Add_rxns,  
                min_threshold, normalize_by_Catoms=normalize_by_Catoms, lb=lb)

            """After confirming that adding the randomly chosen reaction satisfies constraints,
               this reaction is truly appended to the model."""
            M.add_reactions(Add_rxns)
            # Changing possible reaction list for addition and deletion
            possible_adds.remove(add_r)
            possible_dels.append(add_r)

            add_result = ["add_" + add_r] + [str(x) for x in growth_data]
            # Writing to .csv file
            if os.path.exists(csvname):
                with open(csvname,mode='a') as f:
                    f.write(','.join(add_result) + '\n')

            trial = trial + 1
            """If the model reaches to threshold, the knock out experiment is ended."""
            if Avail_C == min_threshold + 1 and not RunToFinal:
                print("The model reach to the threshold...")
                break
        
        cobra.io.write_sbml_model(M,join(model_save_dir,filename + '.xml'))


def long_term_addition_deletion_multi(N):
    expname = secrets.token_hex(4)
    long_term_addition_deletion(model_e_cafba,model_u_e_cafba,Sugar_list,possible_dels,possible_adds,
                                randomwalks,len(Sugar_list)-1,save_folder,expname,
                                RunToFinal = True, normalize_by_Catoms=True, lb=-120.0)
    print("%i process finished"%N)



if __name__ == "__main__":

    homedir = Path(__file__).parent.absolute()
    
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("-m","--model", type=str, default="None", help="an ancestral model for random walks")
    p.add_argument("-u","--universal",type=str,default="None", 
                   help="an universal model (This will be a pool of reactions for adding)")
    p.add_argument("-s","--substrates", type=str, default="None",
                   help="growth substrates (should be given as a list of 'EX_' reactions by txt file)")
    p.add_argument("-r","--randomwalks", type=int, default=10000, help="the number of random walks per trajectory")
    p.add_argument("-n","--Ntrials", type=int, default=1000, help="the number of trials")
    p.add_argument("-t","--threads", type=int, default=1, help="the number of threads")
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
        model_e_cafba = CF.CAFBA_mod_Model(
        cobra.io.read_sbml_model(
            join(homedir,'reference models','iJO1366_cafba.xml')
            )
        )

        model_e = CF.CAFBA_mod_Model(
            cobra.io.read_sbml_model(
                join(homedir,'reference models', 'iJO1366.xml'
                    )
            )
        )
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
        Sugar_list = ['EX_glc__D_e', 'EX_fru_e', 'EX_man_e','EX_fuc__L_e', 
                      'EX_melib_e','EX_gal_e', 'EX_rib__D_e']
    
    else:
        substratefile = args.substrates
        with open(substratefile,mode="r") as f:
            Sugar_list = [x.strip("\n") for x in f.readlines()]


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
    else:
        universalfile = args.universal
        ucafbafile = join(Path(universalfile).parent.absolute(),
                          basename(universalfile).replace(".xml","_cafba.xml"))
        if universalfile.split(".")[-1] == "xml":
            model_u_e = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(universalfile))
            if not os.path.exists(ucafbafile):
                CF.main(universalfile)
            model_u_e_cafba = CF.CAFBA_mod_Model(cobra.io.read_sbml_model(ucafbafile))
        else:
            print("[ERROR] an universal file should be given as .xml file format")
   

    '''Identifying whether each reaction in the model is deletable (non-essential, 
    non-demand,...)'''
    model_e_reactions = {x.id for x in model_e.reactions}
    non_del_reactions = model_e.demand_reactions()
    #non_del_reactions = set(essential_reactions | non_del_reactions)
    non_del_reactions.add('ATPM')
    # possible_dels: deletable reactions in long-term deletion experiment
    possible_dels = set(model_e_reactions - non_del_reactions) & {x.id for x in model_e_cafba.reactions}
    

    print("searching possible reactions for random walks...")
    '''set of the reactions in universal model'''
    total_reactions_id = {x.id for x in model_u_e.reactions}
    sink_rxns = {x.id for x in model_u_e.sinks}
    demand_rxns = {x.id for x in model_u_e.demands}
    exchange_rxns = {x.id for x in model_u_e.exchanges}

    possible_adds = total_reactions_id - set(sink_rxns | demand_rxns | exchange_rxns)
    possible_adds = (possible_adds - model_e_reactions) & {x.id for x in model_u_e_cafba.reactions}
    
    # Here I also excluded transporter reactions in a universal model from this experiment.
    possible_adds = possible_adds - set(model_u_e.transporters_list())
    #print("%i reactions are only found in universal and not in original iJO1366" % len(adding_reactions))

    
    date = str(datetime.datetime.now().date())
    if not os.path.exists(join(os.getcwd(),date)):
        os.mkdir(join(os.getcwd(),date))

    if not os.path.exists(join(os.getcwd(),date,'results')):
        os.mkdir(join(os.getcwd(),date,'results'))

    save_folder = join(os.getcwd(),date,'results')
    model_save_dir = join(Path(save_folder).parent.absolute(),'g_models')
    if not os.path.exists(model_save_dir):
        os.mkdir(model_save_dir)
    #the number of addiiton-deletion trials for making one trajectory
    randomwalks = int(args.randomwalks/2)
    
    #the number of trials of random walks
    Ntrials = args.Ntrials
    
    # # Start the deletion experiments

    with Pool(processes=args.threads) as pool:
        pool.map(long_term_addition_deletion_multi,range(Ntrials))
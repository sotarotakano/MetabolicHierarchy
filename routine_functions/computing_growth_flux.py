__author__ = 'SotaroTakano(sou_tacano@gmail.com)'
__version__ = '1.0.0'

'''Functions for computing growth rate and flux distribution using cobrapy.
All functions are mainly applicable for the simulation under the condition of minimal media (M9)'''

import pandas as pd
import numpy as np
import sys

def get_Catoms_number(met_formula):
    Cpos = met_formula.find('C')
    #If there is no C atoms in this molecule, this normalization is skipped, and None is substituted
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



def computing_flux(EX_C_reactions, model, lb = -20.0, normalize_by_Catoms = False):
    '''compute flux on a given set of carbon sources (EX_C_reactions)'''
    '''return a pandas dataframe of fluxdata in a given model and a give set of Csources'''
    ''' Here, "model" should be CAFBA_mod_model for simplification of the derivating demand reactions'''
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
                M.reactions.get_by_id(i).lower_bound = lb
                growth_rate = M.optimize().objective_value
                if not np.isnan(growth_rate):
                    value = [x.x for x in M.reactions]
                    ID = [x.id for x in M.reactions]
                    fluxdata = pd.Series(value, index = ID)
                    fluxdata.name = i
        if 'flux_C' not in locals():
            flux_C = fluxdata
        else:
            flux_C = pd.concat([flux_C,fluxdata],axis = 1)           
    fluxes_incell = {x.id for x in model.reactions} - model.demand_reactions()
    flux_C  = flux_C.loc[fluxes_incell,:] 
    return flux_C


def computing_flux_pfba(EX_C_reactions, model, lb = -20.0, normalize_by_Catoms = False, internal_reactions = True):
    '''compute flux on a given set of carbon sources (EX_C_reactions)'''
    '''return a pandas dataframe of fluxdata in a given model and a give set of Csources'''
    from cobra.flux_analysis import pfba 
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
                fluxdata = pfba(M).fluxes
                fluxdata.name = i
        if 'flux_C' not in locals():
            flux_C = fluxdata
        else:
            flux_C = pd.concat([flux_C,fluxdata],axis = 1)           
    
    if 'flux_C' not in locals():
        sys.stderr.write("No flux data is generated from a given set of EX_reactions ...\n")
        return None
    else:
        if internal_reactions:
            rxn_incell = {x.id for x in model.reactions} - model.demand_reactions()
            flux_C  = flux_C.loc[rxn_incell,:] 
        return flux_C

def computing_minimal_reaction(EX_C_reactions, model, internal_reactions = True):
    '''compute flux on a given set of carbon sources (EX_C_reactions)'''
    '''return a pandas dataframe of fluxdata in a given model and a give set of Csources'''
    from cobra.flux_analysis import pfba 
    min_rxn = set()
    model.set_minimal_media()
    for i in EX_C_reactions: 
        with model as M:
            if i in M.reactions:
                M.reactions.get_by_id(i).lower_bound = -10.0
                solution = pfba(M)
                growth_rate = solution['Growth']
                if not growth_rate is None and not np.isnan(growth_rate):
                    active_flux = solution.fluxes[abs(solution.fluxes) > 0].index 
                    min_rxn.update(active_flux)      
    if internal_reactions:
        rxn_incell = {x.id for x in model.reactions} - model.demand_reactions()
        min_rxn  =  rxn_incell & min_rxn
    return min_rxn

def jaccard_fluxdata(fluxdata):
    '''compute jaccard distance between every set of fluxdata in a dataframe'''
    Dist = np.array([])
    for i in range(0,fluxdata.shape[0]):
        flux_i = fluxdata.iloc[i,:].to_numpy()
        active_flux_i = {x for x in fluxdata.iloc[i,abs(flux_i) > 1e-06].index}
        for j in range(i+1,fluxdata.shape[0]):
            flux_j = fluxdata.iloc[j,:].to_numpy()
            active_flux_j = {x for x in fluxdata.iloc[j,abs(flux_j) > 1e-06].index}
            j_dist = 1-(len(active_flux_i & active_flux_j)/len(active_flux_i | active_flux_j))
            Dist = np.insert(Dist,len(Dist),j_dist) 
    return Dist

#!/usr/bin/env python
__author__ = 'Sotaro Takano'
__version__ = '1.0.1'

# Updated: Apr,7 2023

''' Extension of CAFBAFY (Jean Villa) to ease the analysis
specifically for the carbon source preference rank (Takano et al., 2023)'''

import sys
import cobra
import os
from .CAFBAFY import * # The original code by Jean Villa
from .computing_growth_flux import *

class CAFBA_mod_Model(CAFBA_Model):
        
    def carbon_metabolites(self):
        '''return the list of carbon metabolites in the model'''       
        exclusion = {'Ca','Co','Cu','Cl','Cd','CrO4'}
        Reactions_EX = set()
    
        Met_C = {x.id for x in self.metabolites if 'C' in x.formula 
                 and not x.formula in exclusion}
    
        #If there is 'EX_' reactions for importing the metabolite of interest, 
        #that metabolite is regarded as Csource.
        for i in Met_C:
            Reactions_i = [x for x in self.metabolites.get_by_id(i).reactions]
            for j in Reactions_i:
                if 'EX_' in j.id:
                    Reactions_EX.add(j.id)
        return Reactions_EX
    
    def carbon_sources(self):
        '''return the list of carbon sources: those of which potentially allow the growth of the model
        when set as a sole carbon source in the environment'''       
        self.set_minimal_media()
        growth_on_Csources = computing_growth(self.carbon_metabolites(), self, 'all')
        C_sources = {x for x in growth_on_Csources[growth_on_Csources > 1e-06].index}
        #print("Possbile carbon sources in the model: %d" % len(C_sources))
        return C_sources


    def essential_reactions(self,EX_reactions):
        '''Screening essential reactions.'''
        '''If single knock-out of a reaction allows no growth on any carbon sources of interest,'''
        '''this reaction is considered as essential.'''  

        from cobra.flux_analysis import single_reaction_deletion
        
        with self as M:
            M.set_minimal_media()
            for k in EX_reactions: 
                if k in M.reactions:
                    M.reactions.get_by_id(k).lower_bound = -10.0
            result = single_reaction_deletion(M)
        essential_rxn = set()
        for i in result[result['growth'] < 1e-06].index:
            essential_rxn.add(list(i)[0])
        print('%i essential reactions in the model' % len(essential_rxn))
        return essential_rxn

    
    def demand_reactions(self):
        '''Screening diffusion, molecular-exchange, demand, and sink reactions.
        All of those are demand reactions in a broad sence especially in our
        randomwalks experiments. Deletion of those reactions are not favourable 
        when we focus on changes in intracellular metabolic network'''
        ex_sink_rxn = {x.id for x in self.reactions if len(x.metabolites) == 1}
        growth_rxn = {x.id for x in self.reactions if x.objective_coefficient == 1}
        trp_rxn = set(self.transporters_list())
        demand_rxn = set(ex_sink_rxn | trp_rxn | growth_rxn)
        return demand_rxn


    def deletable_reactions(self,EX_reactions,additional = set()):
        '''Searching for deletable reactions from a model'''
        '''Diffusion, sink, and essential reactions are eliminated from the deletable reaction pool'''
        '''Additional set of reactions are also excluded if necessary''' 
    
        Reactions = {x.id for x in self.reactions} #All reactions in a given model

        essential = self.essential_reactions(EX_reactions)
        demand = self.demand_reactions()
   
        eliminate_reactions = set(essential | demand | additional)

        KO_reactions = set(Reactions - eliminate_reactions)
        print('%i deletable reactions in the model' % len(KO_reactions))
        return KO_reactions
    
def main(xml_filename,output_file=None,AUTO_Univ=False,AUTO_GLCDe=True):
    x = CAFBA_Model(cobra.io.read_sbml_model(xml_filename))
    x.reassign_compartments()
    x.irreversabilize()
    x.CAFBAfy()
    x.prune_model(AUTO_Univ=AUTO_Univ, AUTO_GLCDe=AUTO_GLCDe)
    x.set_full_media() 
    if output_file==None:
        filename, file_extension = os.path.splitext(xml_filename)
        x.save(filename + '_cafba' + file_extension)
    else: 
        x.save(output_file)

if __name__ == "__main__":
	if len(sys.argv) == 3:
		main(sys.argv[1],sys.argv[2])
	else:
		main(sys.argv[1])
# MetabolicHierarchyEvolution
Python scripts simulating and analyzing evolution of metabolic hierarchy

## Overview	
Briefly, this program generates randomly mutated metabolic models by additions and deletions of metabolic reactions
starting from the user-defined model. The pipeline also tracks the growth capacity and metabolic flux during 
random walks and further extracts evolutionary and metabolic traits.
The primary target of this pipeline is the flexibility of preference rank of available carbon sources, which were
thoroughly studied and documented in Takano et al., 2023 biorxiv (doi:).
All the codes are compatible with cobrapy (Ebrahim et al., 2013 BMC Systems Biol.), 
BiGG models (Norsigian et al., 2020 Nucleic Acids Res.), and CAFBA (Mori et al., 2016 PLoS Comput. Biol.). 

## Requirements
Those packages can be installed by pip or anaconda.

- python (> 3.7)

- cobra

You need to download 

- CAFBAFY.py ([jccvila](https://github.com/jccvila)).

and place it in ./routine_functions. It's available at https://github.com/jccvila/Estrelaetal2021/tree/main/FBA_Simulations.

## Usage
### randomwalks.py
The users can first start from 'randomwalks.py' to generate randomly mutated models. 
The ancestral model can be set by -m option (default is iJO1366) and the universal model by the -u option. 
Universal model is used as a reaction pool for adding the reactions, the default is CUiJO1366.
If you want to generate 100 mutated models by 1000 random walks using 5 threads, run the program as follows.

```bash
python randomwalks.py -n 100 -r 1000 -t 5
```

This creates the mutated models in './[date]/g_models', and the changes in the growth rate of mutated models 
during random walks are saved in './[date]/results'. Those results are used for the following analysis.

### GrowthRank.py
Converting the changes in growth rate to the changes in preference rank during randomwalks.
```bash
python GrowthRank.py -i [result directory]
```
Running this code returns the preference rank changes for a given set of results in './[date]/rank result'.

(You can run this and the following programs without -i option, resulting in processing './test' data (part of results in our study).)

In the original study, one important evolutionary trait is the flexibility of rank flips, 
so the script also returns changes in the preference rank in a given set of substrates
Saved to './[dir]/rankchange_individual_sugars.csv'. 
It also computed the number of rank flip events during each randomwalks experiment. 
Saved to './[dir]/flipcount_pairwise/rankflip_[expid].pkl'

### JaccardPathway.py
Evolutionary changes in metabolic traits are also the primary target of the original study. 
We specifically focused on the dissimilarity in processing pathways across substrates during randomwalks.
You can see them by computing Jaccard distance in used metabolic pathways for all possible pairs of substrates.

```bash
python JaccardPathway.py -r [result directory] -i [interval]
```

This returns the Jaccard distance in the processing metabolic pathways in a given interval in randomwalks. 

### ModifierScreen.py
Screen the key reactions governing the evolvability of the preference rank by comparing the propensity of
rank flips in between the presence and absence of target reactions.
Running the following command generates a summary file of evolutionary modifiers in a given set of results 
Results are summarized into './[dir]/Summary_modifiers_individual.csv'.

```bash
python ModifierScreen.py -i [result directory]
```

### FluxSensitivity.py
The impact of the screened evolutionary modifiers on the metabolic flux is analyzed by 'FluxSensitivity.py'. 
This program computed the flux distribution after deleting or adding the reaction of interest. 
The results are saved into './[dir]/flux_original' and './[dir]/flux_sensitivity_each'

```bash
python FluxSensitivity.py -i [result dir]
```

### ComputeTheta.py
This code computes the flux sensitivity score 'theta' based on the flux distribution files generated by Flux Sensitivity.py.       

```bash
python ComputeTheta.py -i [result dir]
```

### ComputePhi.py
The flux sensitivity on individual metabolite (corresponding to 'phi' in the original paper) is computed by this code.

```bash
python ComputePhi.py -i [result dir]
```

### GetSensitiveMets.py
Screening significantly sensitive metabolites against the mutations on evolutionary modifiers.
The phi scores calculated by ComputePhi.py are used for this statistical screening.

```bash
python GetSensitiveMets.py
```

### RankCarveMeModel.py
You can also compute the preference rank in the prokaryotic models originally generated by CarveMe 
(Machado et al., 2018 Nucleic Acids Res.). 

```bash
python RankCarveMeModel.py
```



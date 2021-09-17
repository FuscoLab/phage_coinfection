# Preventing co-infection: a viral strategy with short-term benefits and long-term drawbacks
Codes for "Preventing co-infection: a viral strategy with short-term benefits and long-term drawbacks" by Michael Hunter and Diana Fusco

- [bioXiv preprint](https://doi.org/10.1101/2021.09.07.458886)

Impotant Note: The code in this repository is currently not well commented (if at all). The scripts will also be a bit of a mess.


## Sections:

- [Abstract](#abstract)
- [Stochastic simulations](#stochastic-simulations)
- [Analysis](#analysis)
- [ODE model](#ode-model)


## Abstract

Viral co-infection occurs when multiple distinct viral particles infect the same host. This can impact viral evolution through intracellular interactions, complementation, reassortment and recombination. In nature many viral species are found to have a wide range of mechanisms to prevent co-infection, which raises the question of how viral evolution is impacted by this strategic choice. Here, we address this question in a model viral system, the ubiquitous bacteriophage and its host bacteria. Using a stochastic model of phage-host interactions in agent-based simulations, we first characterise the behaviour of neutral mutants and find that co-infection decreases the strength of genetic drift. We then quantify how variations in the phage life history parameters affect viral fitness. Importantly, we find that the growth rate (dis)advantage associated with variations in life history parameters can be dramatically different from the competitive (dis)advantage measured in direct-competition simulations. Additionally, we find that co-infection facilitates the fixation of beneficial mutations and the removal of deleterious ones, suggesting that selection is more efficient in co-infecting populations. We also observe, however, that in populations which allow co-infection, a mutant that prevents it displays a substantial competitive advantage over the rest of the population, and will eventually fix even if it displays a much lower growth rate in isolation. Our findings suggest that while preventing co-infection can have a negative impact on the long-term evolution of a viral population, in the short-term it is ultimately a winning strategy, possibly explaining the prevalence of phage capable of preventing co-infection in nature.


## Stochastic simulations

All of the scripts required to run the stochastic simulations used in this work are contained in the `Stochastic_sims` directory.

The scripts can measure the following properties:
1) decay of heterozygosity: `co_het.py` `nonco_het.py`
2) the phage growth rate: `growth_rate.py`
3) the phage fitness in a competative setting: `co_scomp.py` `nonco_scomp.py`
4) the fixation probability of individual mutants: `Wellmixed_co.py` `Wellmixed_nonco.py` `Wellmixed_noncoM_coWT.py` `Wellmixed_noncoWT_coM.py`

In all of the above scripts, `co` and `nonco` indicate whether the population allows or prevents co-infection respectively. In the simulations where the mutant `M` and resident `WT` phage differ in their ability to prevent co-infection, this is also indicated. 

The scripts can be run using the following input parameters:
1) `python (non)co_het.py A B0 beta alpha tau repeat`
2) `python growth_rate.py A B0 beta alpha tau repeat`
3) `python (non)co_scomp.py A B0 beta_M beta_WT alpha_M alpha_WT tau_M tau_WT r0 dr0 repeat`
4) `python Wellmixed_(non)co.py A B0 beta_M beta_WT alpha_M alpha_WT tau_M tau_WT repeat repeat_sub`
5) `python Wellmixed_noncoWT(M)_coM(WT).py A B0 beta_M beta_WT alpha_M alpha_WT tau_M tau_WT repeat repeat_sub`

In all of the above, the symbols have the following meanings:
1) `A`: random number seed (array index for array submissions)
2) `B0`: initial number of bacteria
3) `beta`: phage burst size
4) `alpha`: phage adsorption rate
5) `tau`: phage lysis time
6) `r0`: exponential growth rate of phage at low density
7) `dr0`: error on r0
8) `repeat`: number of simulation repeats 
9) `repeat_sub`: the number of repeats per steady-state simulated

The output of all of the scripts is either printed to the terminal or in a text file, and have the following form:
1) `(non)co_het.py`: `mean heterozygosity`; `fraction of sims with non-zero heterozygosity`
2) `growth_rate.py`: `mean growth rate`, `error on mean`
3) `(non)co_scomp`: `selective (dis)advantage`, `error on mean`
4) `Wellmixed_(non)co.py`: `V_ss`; `dV_ss`; `M_0`; `f_0`; `n_fix`; `P_fix`; `T_fix`
5) `Wellmixed_noncoWT(M)_coM(WT).py`: same as 4)

The results in each row of 4) and 5) represent the average of the `repeat_sub` instances per steady-state. In the above, the symbols have the following meanings:
1) `V_ss`: steady-state number of phage
2) `dV_ss`: error on V_ss
3) `M_0`: number of mutants introduced
4) `f_0`: initial frequency of mutant in free phage population
5) `n_fix`: number of mutant fixation events
6) `P_fix`: frequency of mutant fixation events
7) `T_fix`: mean time to mutant fixation




## Analysis

These scripts can be used to analyse the outputs of the heterozygosity and stochastic simulations described above.

`data_coinfection.mlx` can be run in a folder containing all of the `.txt` files output by `Wellmixed_(non)co.py` and `Wellmixed_noncoWT(M)_coM(WT).py`. It will then save the all of the information in the multiple `.txt` files in one `.mat` file. This file can then be analysed using `Results_coinfection.mlx` which will output the the mean `V_ss`, `f_0`, `P_fix` and their associated errors from all of the simulations carried out.

Similar functions are performed by the `Het_analysis.mlx` and `Het_plot.mlx` scripts for the heterozygosity data. `data_coinfection.mlx` can be run in a folder containing all of the `.txt` files output by `(non)co_het.py`. It will then save the all of the information in the multiple `.txt` files in one `.mat` file. This file can then be analysed using `Het_plot.mlx` which will output the effective population size based on the data and parameters defined in the live script.


## ODE model

Code for ODE

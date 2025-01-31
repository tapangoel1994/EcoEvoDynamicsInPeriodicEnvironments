# *Code for:* Eco-evolutionary dynamics of temperate phages in periodic environments

*By Tapan Goel, 2024*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14782721.svg)](https://doi.org/10.5281/zenodo.14782721)

## General Description

This repository contains all the code need to replicate the figures and analysis show in the paper *Eco-evolutionary dynamics of temperate phages in periodic environments*. The code is written in MATLAB 2022b and MATLAB 2023b. Some functions use parfor loops and therefore need the MATLAB Parallel Computing Toolbox to run. Figure 1 was created in Inkscape 1.3.2 (091e20e, 2023-11-25, custom). Inkscape was also used to modify placement of annotations and legends in figures for clarity.


## Running the code

There is no existing data being used in this paper. All the results are based on solving the model setup in the methods section of manuscript. You will need MATLAB 2022b or MATLAB 2023b to run all the code (works in both versions). You will need to have the MATLAB Parallel Computing Toolbox to run the functions "PopulationSteadyStateFunction", "InvasionDynamics", "Stochastic_QL_Stochastic_QV_highdeathrate" and "Stochastic_QL_Stochastic_QV_highdeathrate". I have hardcoded the use of 64 CPU nodes for parallelizing the computation involved in generating those figures. You can change that number for running the code on your machine, by changing the NumNodes variable in the script that generates the figure you want to reproduce. To generate a particular figure, set ./Code as your working directory and run the appropriate script to generate the corresponding figure.

If you're running the code on mac/linux instead of windows, you might have to change the file separator ("/" for mac/linux and "\\" for windows). You will need to do this for all the figure generating scripts and the functions in the utils folder.

## Notes for code reviewer

I suggest reviewing the code files in the following order:
   1. Read the methods section, and Tables S1 and S2 of the manuscript.
   2. Match Equation 7 in the manuscript with the equations in Code\utils\ODE_RSEILV_2species.m
   3. Match the parameter values in Tables S1 and S2 in the manuscript with the numbers in the file Code\lib\fixedparameters.m
   4. Read through the code in Code\utils\PopulationSteadyStateFunction.m
   5. Read through the code in Code\utils\InvasionDynamics.m

Now, you should be able to generate all the figures in the paper (except figure 1) by simply running the corresponding script file. Make sure you are in the ./Code directory while running those scripts.

## Folder content description

### ./Doc

Contains the most up-to-date version fo the manuscript.

### ./Data

Contains .mat files generated by running the scripts Figure4.m, Figure5.m and Figure6.m.

- **"SteadyState_CyclePeriod=\<CyclePeriod\>,d=\<celldeathrate\>,q_L=\<q_L\>,q_V=\<q_V\>.mat"**: contains steady state viral densities and number of cycles required to reach steady state for a range of integration probabilities and induction rates, given a growth cycle duration of *Cycle Period*, cell death rate *d*, lysogen passaging fraction *q_L* and virion passaging fraction *q_V*.

- **"Invasion_CyclePeriod=\<CyclePeriod\>,d=\<celldeathrate\>,Gamma=\<Gamma\>,q_L=\<q_L\>,q_V=\<q_V\>.mat"**: contains viral densities at invasion, matrix corresponding to the pairwise invasion plots and number of cycles required to ascertain invasion success or failure, given a growth cycle duration of *Cycle Period*, cell death rate *d*, gamma *gamma*, lysogen passaging fraction *q_L* and virion passaging fraction *q_V*.

- **Stochastic_QL_Stochastic_QV_highdeathrate.mat :** contains number of cycles each viral strategy survives in each stochastic realization with different ranges of q_L and q_V for the high cell death rate case.

- - **Stochastic_QL_Stochastic_QV_lowdeathrate.mat :** contains number of cycles each viral strategy survives in each stochastic realization with different ranges of q_L and q_V for the low cell death rate case.

### ./RevisedFigures

Output folder for figures. Contains subfolder for the main text figures, SI figures associated with the main text figures and SI figures associated with the high cell death rate case.

### ./Code

Contains all code associated with the manuscript.

- **Code/lib/colorpalette.m :** Contains font and color scheme settings used in figures throughout the manuscript. This script is run at the beginning of each figure generating script.

- **Code/lib/fixedparameters_lowdeathrate.m :** Contains all parameters that are help constant across different simulations associated with the manuscript for the low cell death rate case. This script is run at the beginning of each figure generating script.

- **Code/lib/fixedparameters_highdeathrate.m :** Contains all parameters that are help constant across different simulations associated with the manuscript for the high cell death rate case. This script is run at the beginning of each figure generating script.

- **Code/lib/VariableMaps.m:** Contains parameters that map variable names to vector indices

- **Code/utils/ODE_RSEILV_2species.m :** Returns derivative of system dynamics at time t given life history parameters and system state at time t. Computes *Equation (7)* presented in the main manuscript.

- **Code/utils/PopulationSteadyStateFunction.m :** Function computes the steady state of the virus host system for given values of life history traits, initial conditions, growth cycle duration and filtration conditions. Output gets stored in "../Data/SteadyState_CyclePeriod=\<CyclePeriod\>,d=\<celldeathrate\>,q_L=\<q_L\>,q_V=\<q_V\>.mat". *Note:* Function uses a parfor loop and takes the number of nodes over which to parallelize as input.

- **Code/utils/InvasionDynamics.m :** Function computes pairwise invasions between resident-mutant pairs with given life history traits, initial conditions, growth cycle duration and filtration conditions. Output gets stored in "../Data/Invasion_CyclePeriod=\<CyclePeriod\>,d=\<celldeathrate\>,Gamma=\<Gamma\>,q_L=\<q_L\>,q_V=\<q_V\>.mat". *Note:* Function uses a parfor loop and takes the number of nodes over which to parallelize as input.

- **Code/utils/FindESSStrategy.m :** Function to calculate the ESS strategy (if one exists) given a pairwise invasibility plot and the phenotype thats changing.

- **Code/utils/StochasticRealization.m :** Function to calculate the number of cycles until which a virus with a given strategy survives in a stochastic environment

- **Code/DataGenerate_highdeathrate.m :** Script loops over filtration conditions and cycleperiods to generate steady state densities for a range of viral strategies for the high cell death rate case. Generates the "../Data/SteadyState_CyclePeriod=\<CyclePeriod\>,d=\<celldeathrate\>,q_L=\<q_L\>,q_V=\<q_V\>.mat" files.

- **Code/DataGenerate_lowdeathrate.m :** Script loops over filtration conditions and cycleperiods to generate steady state densities for a range of viral strategies for the low cell death rate case. Generates the "../Data/SteadyState_CyclePeriod=\<CyclePeriod\>,d=\<celldeathrate\>,q_L=\<q_L\>,q_V=\<q_V\>.mat" files.

- **Code/Stochastic_QL_Stochastic_QV_highdeathrate.m :** Script loops over different ranges of stochastic qV and stochastic qL to generate a table of number of cycles that different strategies survive in several iterations of the numerical experiment. Generates the file "../Data/Stochastic_QL_Stochastic_QV_highdeathrate.mat".

- **Code/Stochastic_QL_Stochastic_QV_lowdeathrate.m :** Script loops over different ranges of stochastic qV and stochastic qL to generate a table of number of cycles that different strategies survive in several iterations of the numerical experiment. Generates the file "../Data/Stochastic_QL_Stochastic_QV_lowdeathrate.mat".

- **Figure2.m :** generates figure 2.

- **Figure3.m :** generates figure 3.

- **Figure3_SI_PopDynamics.m :** generates figures S1-S3 and S6-S8.

- **Figure3_SI_FinalStateHistogram_lowdeathrate.m :** generates figure S4.

- **Figure3_SI_FinalStateHistogram_highdeathrate.m :** generates figure S9.

- **Figure4_lowdeathrate.m :** generates figure 4.

- **Figure4_highdeathrate.m :** generates figure S10.

- **Figure5_lowdeathrate.m :** generates figure 5.

- **Figure5_highdeathrate.m :** generates figure S11.

- **Figure6_lowdeathrate.m :** generates figure 6.

- **Figure6_highdeathrate.m :** generates figure S13.

- **CycleToCycleInvasionDynamics.m :** generates Figure S12.


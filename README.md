# *Code for:* Eco-evolutionary dynamics of temperate phages in periodic environments

*By Tapan Goel, 2024*

## General Description

This repository contains all the code need to replicate the figures and analysis show in the paper *Eco-evolutionary dynamics of temperate phages in periodic environments*. The code is written in MATLAB 2022b and MATLAB 2023b. Some functions use parfor loops and therefore need the MATLAB Parallel Computing Toolbox to run. Figure 1 was created in Inkscape 1.3.2 (091e20e, 2023-11-25, custom).

## Folder content description

### ./Data

Contains .mat files generated by running the scripts Figure4.m, Figure5.m and Figure6.m.

- **"SteadyState_CyclePeriod=\<CyclePeriod\>,S0=\<InitialHostDensity\>,V0=\<InitialViralDensity\>,q_L=\<q_L\>,q_V=\<q_V\>.mat"**: contains steady state viral densities and number of cycles required to reach steady state for a range of integration probabilities and induction rates, given a growth cycle duration of *Cycle Period*, initial host density *S0*, initial viral density *V0*, lysogen passaging fraction *q_L* and virion passaging fraction *q_V*.

- **"Invasion_CyclePeriod=\<CyclePeriod\>,S0=\<InitialHostDensity\>,V0=\<InitialViralDensityofResident\>,q_L=\<q_L\>,q_V=\<q_V\>.mat"**: contains viral densities at invasion, matrix corresponding to the pairwise invasion plots and number of cycles required to ascertain invasion success or failure, given a growth cycle duration of *Cycle Period*, initial host density *S0*, initial viral density *V0*, lysogen passaging fraction *q_L* and virion passaging fraction *q_V*.

### ./Figures

Output folder for figures generated by figure producing scripts.

### ./Code

Contains all code associated with the manuscript.

- **Code/lib/colorpalette.m :** Contains font and colorscheme settings used in figures throughout the manuscript. This script is run at the beginning of each figure generating script.

- **Code/lib/fixedparameters.m :** Contains all parameters that are help constant across different simulations associated with the manuscript. This script is run at the beginning of each figure generating script.

- **Code/utils/ODE_RSEILV_2species.m :** Returns derivative of system dynamics at time t given life history parameters and system state at time t. Computes *Equation (7)* presented in the main manuscript.

- **Code/utils/PopulationSteadyStateFunction.m :** Function computes the steady state of the virus host system for given values of life history traits, initial conditions, growth cycle duration and filtration conditions. Output gets stored in "../Data/SteadyState_CyclePeriod=\<CyclePeriod\>,S0=\<InitialHostDensity\>,V0=\<InitialViralDensity\>,q_L=\<q_L\>,q_V=\<q_V\>.mat". *Note:* Function uses a parfor loop and takes the number of nodes over which to parallelize as input.

- **Code/utils/InvasionDynamics.m :** Function computes pairwise invasions between resident-mutant pairs with given life history tratis, initial conditions, growth cycle duration and filtration conditions. Output gets stored in "../Data/Invasion_CyclePeriod=\<CyclePeriod\>,S0=\<InitialHostDensity\>,V0=\<InitialViralDensity\>,q_L=\<q_L\>,q_V=\<q_V\>.mat". *Note:* Function uses a parfor loop and takes the number of nodes over which to parallelize as input.

- **Code/Figure2.m :** Generates figure 2.

- **Code/Figure3.m :** Generates figure 3 in the manuscript. The diagrams in panels B,D and F in figure 3 were added to the plots generated by Figure3.m later using Inkscape.

- **Code/Figure4.m :** Generates figure 4 in the manuscript. Uses the PopulationSteadyStateFunction and InvasionDynamics functions and therefore needs the Parallel Computing ToolBox.

- **Code/Figure5.m :** Generates figure 5 in the manuscript. Uses the PopulationSteadyStateFunction and InvasionDynamics functions and therefore needs the Parallel Computing ToolBox.

- **Code/Figure6.m :** Generates figure 6 in the manuscript. Uses the PopulationSteadyStateFunction and InvasionDynamics functions and therefore needs the Parallel Computing ToolBox.

- **Code/Figure7.m :** Generates figure 7 in the manuscript.

- **Code/FigureS1_2_3.m :** Generates SI figures 1,2 and 3.

- **Code/FigureS4.m :** Generates SI figure 4.

- **Code/FigureS5.m :** Generates SI figure 5.

- **Code/FigureS6.m :** Generates SI figure 6.

- **Code/FigureS7.m :** Generates SI figure 7.

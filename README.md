# WWT microcolony sizes
This repository contains the codes used to run the scripts for the paper **"Modelling the size distribution of activated sludge microcolonies"**, by Niccolò Totis, An-Sophie Christiaens, Ilse Smets and Steffen Waldherr.

## Getting Started

The script **Main.m** can be used to run all simulations displayed in the paper. The script **Main_lowerSRT_onlyPlots.m** can be used to plot the last figure of model interrogation where lower SRT values are used.

## Simulation options

  - **parSetup** defines the set of kinetic parameters used in the ASM1seq:
    - 'baseline' refers to the parameter set reported in the original paper presenting the ASM1 model [Henze, 2000]
    - 'EEMorris' can be used to run the sensitivity analysis via the method of Morris
    - 'paramSweep' performs ASM1seq parameter sweep procedure
    - 'chosenSetASM1seq' is the identified parameter set, later used to simulate the whole hybrid model (ASM1seq + PBE)
  - **IC_at** defines the initial condition:
    - 'data' implies that the inital state of the model is assigned to that inferred from the steady-state experimental data collected from the system. This option needs to be used during the identification of the ASM1seq model, which, as the Results section of the paper describes, consists of a sensitivity analysis via the method of Morris and a parameter sweep procedure.  
    - 'SS_computed' is used to reproduce the steady-state scenario, where the model is initialised to the (computed) steady state that is most coherent with the data;
    - 'lower_biomass' is used to simulate the hypothetical startup scenario, described in the paper, aimed to assess the model robustness. 

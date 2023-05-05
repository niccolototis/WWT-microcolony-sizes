 function P = GeneralSettings(P)
 %%%%%%%
 % NB IMPO do not delete global variables even if not used otherwise they
 %are not stored in structure P
 %%%%%%%
global  implementGrowth washOutParticles checkECDF tspanCut T_sim_PBE T_sim saveParSweepResToFile debugParallel particleInflux modelStructure IC_at T_SS whatSweep sweep_lambda parSetup identify_lambda biom_lower_than_SS
global  influxActive lambdas cycleModified cycleModStr biomFract0 modelWaterSludgeSegreg computeSS_for_later_IC instantSludgeRemoval variableInGrams kineticParams implementBreakage modelEPS modelI
global runPBE plotSimulation plotStackedSimul plotOnlyEstimatedDistr plotOnlyFirstEndSlices 
global parametersEffectOnMSD noTrials plotLambdaProfile

%% - What parameter set: activate one of the 4 following options
% % parSetup = 'baseline';
% % parSetup = 'EEMorris';
% % parSetup = 'paramSweep';
% % parSetup = 'chosenSetASM1seq';

if ~cycleModified
    cycleModStr = '';
end

if strcmp(IC_at,'data') || strcmp(parSetup,'paramSweep') || strcmp(parSetup,'EEMorris')
    runPBE = false;
    plotSimulation = false; 
    plotStackedSimul = false; 
    plotOnlyEstimatedDistr = false;
    plotOnlyFirstEndSlices = false;
    plotLambdaProfile = false;
    parametersEffectOnMSD = false; noTrials = 0;
end

if P.plotPhases
     parSetup = 'chosenSetASM1seq';
     IC_at = 'SS_computed';
     P.runASM1seq = false;
     P.T_sim=50;
end

if P.plotLambdaProfile
     parSetup = 'chosenSetASM1seq';
     IC_at = 'SS_computed';
     P.runASM1seq = false;
end

switch parSetup
 case 'chosenSetASM1seq'
    P.nChosenParSet_ASM1seq = [];
 case 'EEMorris'
    IC_at = 'data';
    P.plot='both_HET_AUT';
    
    % - Select the type of figure
    P.figure = 'surface';
    
    disp('For the Morris Method I am starting from the baseline parameter values!')
 case 'paramSweep'
     P.numberSweepBatches_ASM1seq = ceil(P.nSweeps/P.batchSizeSweep);
    IC_at = 'data';
    % - Select one of the following options for how param combinations are
    % generated
    whatSweep = 'LHS';
    % whatSweep = 'Ladder';
    % whatSweep = 'meshgrid';
    
    debugParallel = true;
    
    % - Select the type of figure
    P.figure = 'flat';% 'flat'=normal plot of biomass
    % P.figure = 'surface';
    
    % - Select what trajectory is used for the comparison
    % P.plot = 'AUT'; % 
    % P.plot = 'HET'; %
    P.plot = 'both_HET_AUT'; % 
end

%% - choose how the initial conditions are defined 
switch IC_at
 case 'data'
    % in this case I am running parameter sweep to indentify a
    % parameter set that brings the final state to a state that is
    % considered realistic according to the experimentalists. This
    % final (steady) state is then used as the intial condition for
    % subsequent runs, so that the final dynamics are flat (SS is maintained)
 case 'SS_computed'
    % This is the initial condition I use to reproduce the data, and so
    % to estimate the best lambda value
    % Note that only if T_sim = 50 I have actually reached the steady
    % state
 case 'startupScenario'
    % After lambda has been estimated, I initialize the system to a
    % lower biomass and check if the sustem converges to the same
    % steady state. Here I specify the fraction of biomass that is used
    biomFract0 = 1/100; 
end


% - Others
modelStructure = 'actualSystem'; % 'ABMcurtis'
kineticParams = 'actualSystem';%'ASM1_original'; 
variableInGrams = false;
modelWaterSludgeSegreg = 'segregated'; %'stirred';
instantSludgeRemoval = false;
implementGrowth = true; 
implementBreakage = true; 
% --
influxActive = true;
particleInflux = true;
washOutParticles = true;

checkECDF = true; % see in kdestimation

namesGlobal = who('global'); %# A cell array of variable names
namesLocal = who(); 
namesGLocal = intersect(namesGlobal,namesLocal); 
for iVar = 1:numel(namesGLocal)
 P.(namesGLocal{iVar}) = eval(namesGLocal{iVar}); % [EDITED]
end
return

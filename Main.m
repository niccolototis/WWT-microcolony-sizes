clear all
close all
global T_sim parSetup lambdasToSweep figsVisible nDecimalsToIncludeNDFsampling
global timeLastCycleStart repeat234 allChecks step indSwitchTimes phaseCtrl 
global plotOnlyFirstEndSlices runNewSweeps sweep_lambda plotOnlyEstimatedDistr
global labLambda noTrials iT  parametersEffectOnMSD improvStr plotWhat
global cycleModStr howDistancesComputed saveParSweepResToFile  
global distrib_at_startupScenario folderOutput IC_at speciesPBE spNames

% To only plot the impact of reduced SRT (last figure on the paper), check
% Main_final_lowerSRT_onlyPlots

addpath(genpath('./')) % do not delete this, it is needed to access additional scripts 

%% Current main settings 
% P is a structure containing all options
P.runASM1seq = true; P.plotASM1seq = true; 
runPBE = true; 
plotSimulation = true; 
plotStackedSimul = true; 
P.plotPhases = false; 
plotOnlyEstimatedDistr = false;
plotOnlyFirstEndSlices = true;
P.plotLambdaProfile = false;

% Options referring to the last set of analyses described in the paper,
% exploring how changes in controllable parameters would impact the steady
% state MSD in the described WWT plant
parametersEffectOnMSD = false; noTrials = 0;


%% Whole bag model options
P.tryCatch = true;

% Choose the aSpPBEecific setup for the parameters used
parSetup = 'chosenSetASM1seq'; 
% 'baseline' uses the baseline parameter set of ASM1 model in Henze, 2000 
% 'EEMorris' can be used to run the sensitivity analysis via the method of Morris
% 'paramSweep' performs a sweep in parameters 
% 'chosenSetASM1seq' is the set of parameters eventually chosen

% Choose option for the initial condition
IC_at = 'SS_computed'; 
% 'data' is used when running the sensitivity analysis and the parameter sweep
% 'SS_computed' initializes ASM1seq to the final state obtained after the parameter sweep
% 'startupScenario' simulates the startup scenario

% Length of the simulation
T_sim = 50; % days. It should be a multiple of 0.5 = operation cycle length 

% Options for parameter sweep of ASM1seq model
runNewSweeps = false; P.nSweeps = 5;P.batchSizeSweep = 50; P.nBestKparSweeps = 26;
saveParSweepResToFile = true;
if strcmp(IC_at,'data') || strcmp(parSetup,'paramSweep') || strcmp(parSetup,'EEMorris')
    parametersEffectOnMSD = false;
end

%% PBE options
aSpPBE = 'AUT'; % PBE simulates the dynamics of the NDF for just autotrophs
figsVisible = 'on';
computeBestLambda = false;

% plotting
if plotOnlyEstimatedDistr || plotOnlyFirstEndSlices
     plotSimulation = true;
     plotStackedSimul = true;
end

% Options for the comparison of computed vs experimental MSD
plotWhat = 'N'; % Stands for the NDF number density function 
plotWhatprofile.stat = 'CVM';
plotWhatprofile.weight = 'equal';
distrib_at_startupScenario = 'impulse'; %{'lognormal','normal','impulse'}

% Sweeping lambda and finding its optimal value
sweep_lambda = false;
lambdasToSweep = [0.0001 0.001 0.01 0.1 0.5 0.8 [1:0.5:9.5] [10:20] [25:5:120]];
nDecimalsToIncludeNDFsampling = 3;
lambdaVs = 'sum_dist_all'; %{'sum_dist_all','dist_tend'};
statsToCompute = {'KS','AD','CVM'};
logDistance = [false true true];
howDistancesComputed = 'tptsAllin'; 
% With 'tptsAllin' distances between 2 different simulations are computed by 
% summing up distances of all the slices of the simulations
% With 'tptsApart' distances are computed separately between each pair of
% slices and then they are averaged
optimal_lambda = 52;
lambdasToplot = optimal_lambda;

% I run the for loop only if I am solving the PBE for every bestSweep
% parameter set identified
switch parSetup
 case 'chosenSetASM1seq'
     P.nBestKparSweeps = 1; 
     P.filename_pSweepASM1seq_results_toPlot = './outputs/parameter_sweep_ASM1seq/results_paper.mat';
 case 'baseline'
    P.nBestKparSweeps = 1; 
end

P = paramsToModify(P);
% These improvements are changes of the baseline operation parameters
% or kinetic parameters done as the last experiment for the paper, to
% see if it is possible to improve the particle size distribution (i.e.
% shifting it to the right)

for iT = 0:noTrials 
    disp(iT)
     if parametersEffectOnMSD 
        improvStr = ['_impro_' char(datetime('today')) '_v' num2str(iT)];
     else
         improvStr = '';
         noTrials = 0;
     end
     P = GeneralSettings(P); 
     P = defineIndexes(P);
     P = setOperationParams(P);

     %% set and run whole bag model
     if runPBE || P.runASM1seq || P.plotASM1seq
         phaseCtrl = 1;timeLastCycleStart = 0;repeat234 = 1; allChecks = 1; step = 0; indSwitchTimes = 1;% initializing control
         %[P,t_all_ASM1seq,X_c_tall,S_c_tall,V_tall] = runASM1seqModel_ASM1_CURRENT(P);

         [P,t_all_ASM1seq,X_c_tall,S_c_tall,V_tall] = runASM1seq(P);
         % The output of runASM1seqModel_ASM1 is concentrations. Here I
         % revert to grams to compare with the PBE

         % Xgrams_ASM1seq was only actually used in the plotBiom script,
         % that was meant to be compared to the first moments
         Xgrams_ASM1seq = X_c_tall .* V_tall; 

         %% set and run PBE
         % In the following M_N is a 2D matrix, storing the value of n given the two
         % independent coordinate values x1 and x2
         if runPBE

             P = PBEsettings(P);
             plotsVisible = false;
             [aSt,P] = PBEParams(P,plotsVisible); % kinetic parameters, initial conditions, inlet concentrations, PBEparameters (mu sigma) 
             [M_Nti_0_full,P] = defineIC(aSt,P); % Each column of M_Nti_0_full contains [N_all_aChar, x1, w1, Ntot]
             [Ker_B,Ker_V] = defineBrekageKernel(aSt); % Define breakage kernels for all sizes
             
             for aSpPBE=speciesPBE
                sp=char(spNames{aSpPBE}); 
                 % This mu is the connection between ASM1seq and PBE
                 mu_ASM1seq_tall = compute_mu(S_c_tall,sp,P); 
                 if sweep_lambda
                     if debugParallel
                         first = 1; last = length(P.lambdasToSweep);
                         P.figure = 'flat';
                         results = collectResultOnePartitionPBE(first,last,aSt,P,M_Nti_0_full,t_all_ASM1seq,mu_ASM1seq_tall,Ker_B,Ker_V,sp,[]);
                     else
                        parallelRuns_PBE(aSt,P,M_Nti_0_full,t_all_ASM1seq,mu_ASM1seq_tall,Ker_B,Ker_V,sp) 
                     end
                 else
                     P.lambda = optimal_lambda;
                     aSt = runPBEforOneParam(aSt,P,M_Nti_0_full,t_all_ASM1seq,mu_ASM1seq_tall,Ker_B,Ker_V,sp);
                 end
             end             
         end
     end

    %% Compute p-values for simulations produced for different lambdas
    if computeBestLambda
        computeBestLambda(P)
    end

    %% Plot the pval vs lambda profile
    if P.plotLambdaProfile
        plotLambdaProfile(statsToCompute,plotWhatprofile,lambdaVs,howDistancesComputed,divideDistByCritValue,logDistance)
    end

    %%   Plots 
    if plotSimulation
        close all
        switch IC_at
            case 'SS_computed'

                % - plot for each lambda
                for i=1:length(lambdasToplot)
                    Lambda_i=lambdasToplot(i);
                    labLambda=num2str(Lambda_i);
                    fprintf('##### Plotting n=%g for lambda %g ##### \n ',i,Lambda_i)

                    % - load PBE data
                    switch P.parSetup
                        case 'chosenSetASM1seq'
                            load(strcat([P.folderOutput,'P_allLambdas/P_lambda_',labLambda,improvStr,cycleModStr,'.mat']),'P'); 
                            load(strcat([P.folderOutput,'aSt_allLambdas/aSt_lambda_',labLambda,'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells), improvStr, cycleModStr,'.mat']),'aSt');
                    end
                    if ~isfield(P,'initialNDF') && isfield(P,'initialPdf')
                        P.initialNDF=P.initialPdf;
                    end

                    % - load compPD, used to compute the quartiles
                    foldPD=[folderOutput 'S_compPD/']; labLambda=strrep(labLambda,'.','p');
                    ccc=load([foldPD 'S_compPD_lambda_' labLambda '.mat'],'compPD'); 
                    compPD=ccc.compPD;

                    % - prepareToPlot
                    aSt=prepareToPlot(aSt,aSt,compPD);
                    if plotStackedSimul
                        % - plot
                        [joyFigs]=prepareJoyPlot(P,aSt);
                        % - save
                        foldFig='./outputs/figures/joyPlot/';
                        if strcmp(P.parSetup,'bestSweep_ASM1seq_to_PBE')
                            foldFig=[foldFig,'alternative_',vv,'/'];
                            if ~exist(foldFig,'dir')
                                mkdir(foldFig)
                            end
                        end
                        saveFigs(joyFigs,foldFig)
                    end
                end
            case 'lower_biom'
                labLambda=num2str(optimal_lambda);

                % load structures containing options
                load(strcat([P.folderOutput,'P_allLambdas/P_', distrib_at_lower_biom ,'_lambda_',labLambda,improvStr,cycleModStr,'.mat']),'P'); 
                P=PBEsettings(P);
                plotsVisible=false;
                [aSt_options,P]=PBEParams(P,plotsVisible); 

                % - load compPD, used to compute the quartiles
                foldPD=[folderOutput 'S_compPD/']; 
                ccc=load([foldPD 'S_compPD_lambda_' labLambda '.mat'],'compPD'); 
                compPD=ccc.compPD;

                % - plot for each lambda
                fprintf('##### Plotting n=%g for lambda %g ##### \n ',i,optimal_lambda)                    

                % - load PBE data low_biom
                load(strcat([P.folderOutput,'aSt_allLambdas/aSt_', distrib_at_lower_biom ,'_lambda_',labLambda,'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells),'.mat']),'aSt');
                % - load compPD, used to compute the quartiles
                foldPD=[folderOutput 'S_compPD/']; 
                ccc=load([foldPD 'S_compPD_lambda_' labLambda '.mat'],'compPD'); 
                compPD=ccc.compPD;
                % - prepareToPlot
                aSt=prepareToPlot(aSt,aSt_options,compPD);

                % - prepareToPlot
                    aSt=prepareToPlot(aSt,aSt_options,compPD);

                % - load PBE data SS_computed for the same optimal_lambda for joyFigs.fig2 (CFR endpoints)
                folderOutput_SS=strrep(P.folderOutput,'lowBiom','SS');
                aSt_SS_computed=load(strcat([folderOutput_SS,'aSt_allLambdas/aSt_lambda_',labLambda,'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells),cycleModStr,'.mat']),'aSt');
                aSt_SS=aSt_SS_computed.aSt;
                % - load compPD, used to compute the quartiles
                foldPD=[folderOutput_SS 'S_compPD/']; 
                ccc=load([foldPD 'S_compPD_lambda_' labLambda '.mat'],'compPD'); 
                compPD_SS=ccc.compPD;
                % - prepareToPlot
                aSt_SS=prepareToPlot(aSt_SS,aSt_options,compPD_SS);

                % - plot
                if plotStackedSimul
                    joyFigs=prepareJoyPlot(P,aSt,aSt_SS);
                    % - save
                    foldFig='./outputs/figures/joyPlot/';
                    saveFigs(joyFigs,foldFig)
                end
        end
    end
end

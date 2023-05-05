clear all
close all
global T_sim modelStructure nCells checkECDF parSetup lambdasToSweep speciesPBE spNames nTpoints_all nCellsGh nChar nTimeFragms tspans_ASM1seq N_tot_0_all plotPVF plotForPresentation figsVisible
global timeLastCycleStart repeat234 allChecks step indSwitchTimes phaseCtrl plotOnlyFirstEndSlices runNewSweeps thisTpoint steP broken tFragm_lengths sweep_lambda plotOnlyEstimatedDistr
global cellDensity x2_ind_toPlot noTrials iT sampleImprovements parametersEffectOnMSD howDistancesComputed improvStr cycleModified cycleModStr plotCfrEndSlice plotForPresent saveParSweepResToFile plotWhat distrib_at_startupScenario folderOutput initialPdf M_Nti_0_pdf out_vAV IC_at identify_lambda debugParallel nDecimalsToIncludeNDFsampling

% Minimal script to plot the SS scenario sweeping SRT (trying to improve
% the curve, i.e. swithcjing it to the right)

nDecimalsToIncludeNDFsampling = 3;plotForPresentation=true;
P.forVSC = false;
if P.forVSC
 addpath('/apps/leuven/skylake/2018a/software/CPLEX/12.10.0-intel-2018a-Python-3.6.4/cplex/matlab/x86-64_linux')
end
addpath('./2DPBE_src') % do not delete this, it is needed to access additional scripts 

%% Current main settings 
P.runASM1seq = false; P.plotASM1seq=true; computeBestLambda=false; plotSimulation=true;
P.plotSubs = false;

%
parSetup = 'chosenSetASM1seq'; % {'EEMorris' 'chosenSetASM1seq' 'baseline' 'paramSweep'}
parametersEffectOnMSD = true; noTrials=10; cycleModified=false; sampleImprovements=false;

%% Whole bag model options
figs_ASM1seq = [];
P.tryCatch = true;
P.plotPhases = false; % Plotting only the first cycle
IC_at = 'SS_computed'; %{'data','SS_computed','startupScenario'};
T_sim = 50; % days. It should be a multiple of 0.5 = operation cycle length 

% Parameter sweep of ASM1seq model
runNewSweeps = false; 

%% PBE options
figsVisible = 'on';
% plotting
plotStackedSimul = true; joyFigs=[]; 
otherPlots = false;
plotOnlyEstimatedDistr = false;
plotOnlyFirstEndSlices = false;
P.plotLambdaProfile = false;
plotForPresent = false;
if plotOnlyEstimatedDistr || plotOnlyFirstEndSlices
 plotSimulation = true;
 plotStackedSimul = true;
end

%
plotWhat = 'N'; %{'N','pdf','pvf'};
plotWhatprofile.stat = 'CVM';
plotWhatprofile.weight = 'equal';
distrib_at_startupScenario = 'impulse'; %{'lognormal','normal','impulse'}

% sweeping lambda and finding its optimal value
sweep_lambda = false;
lambdaVs = 'sum_dist_all'; %{'sum_dist_all','dist_tend'};
divideDistByCritValue = false;
statsToCompute = {'KS','AD','CVM'};
logDistance = [false true true];
howDistancesComputed = 'tptsAllin'; 
optimal_lambda = 52;
lambdasToplot = optimal_lambda;
lambdasToSweep = optimal_lambda;

% I run the for loop only if I am solving the PBE for every bestSweep
% parameter set identified
switch parSetup
 case 'chosenSetASM1seq'
    P.nBestKparSweeps = 1; 
    P.filename_pSweepASM1seq_results_toPlot = './outputs/parameter_sweep_ASM1seq/results_paper.mat';
end

P = paramsToModify(P);

start_iT = 0; P.allObj=[]; P.allObjNames=[];

% I start from 0, for iT ==0 I plot the results with the actual SRT
for iT = start_iT:noTrials 
 
 if iT>0 
    parametersEffectOnMSD = true;
    % This is the date when I produced those results
    improvStr = ['_impro_09-Dec-2022_v' num2str(iT)];
 else
    parametersEffectOnMSD = false;
    improvStr = '';
 end
 
 P.iSw = 1;
 P = GeneralSettings(P);  
 P = defineIndexes(P);
 P = setOperationParams(P);
 
 % Store names for legend
 P.allObjNames = [P.allObjNames {['$\mathbf{',num2str(P.SRT),'} \, \mathrm{\mathbf{d;}}$']}];
 
 plotToSampleColors(start_iT)       
 [P,figs_ASM1seq] = plot_ASM1seq_lowSRT(P,figs_ASM1seq);
 
 
 %% NB this is necessary to plot the correct statistics for each iT
 if computeBestLambda
    computeBestLambda_fn(P,'for') % {'for','parfor'}
 end
 
 %% Plot simulation
 if plotSimulation
    switch IC_at
        case 'SS_computed'
            aSt_SS_computed = [];
            if strcmp(P.parSetup,'chosenSetASM1seq') && exist([folderOutput,'P_allLambdas/P_lambda_0.0001' improvStr cycleModStr '.mat'],'file')
                
                % load structures containing options
                P_opt = load(strcat([folderOutput,'P_allLambdas/P_lambda_0.0001',improvStr cycleModStr,'.mat']),'P');
                P_opt = PBEsettings(P_opt.P);
                plotsVisible = false;
                [aSt_options,P_opt] = PBEParams(P_opt,plotsVisible); 
            end

            % plot for each lambda
            for i = 1:length(lambdasToplot)
                Lambda_i = lambdasToplot(i);
                labLambda = num2str(Lambda_i);
                fprintf('##### Plotting n = %g for lambda %g ##### \n ',i,Lambda_i)

                % load PBE data
                switch P_opt.parSetup
                    case 'chosenSetASM1seq'
                        oneP = load(strcat([P_opt.folderOutput,'P_allLambdas/P_lambda_',labLambda,improvStr,cycleModStr,'.mat']),'P'); 
                        oneSt = load(strcat([P_opt.folderOutput,'aSt_allLambdas/aSt_lambda_',labLambda,'_T_sim',num2str(P_opt.T_sim),'_nCells',num2str(P_opt.nCells), improvStr, cycleModStr,'.mat']),'aSt');
                end
                if ~isfield(oneP.P,'initialNDF') && isfield(oneP.P,'initialPdf')
                    oneP.P.initialNDF = oneP.P.initialPdf;
                end
                oneP.P.allObj = P.allObj; oneP.P.allObjNames=P.allObjNames; 

                % load compPD, used to compute the quartiles
                foldPD = [folderOutput 'S_compPD/']; labLambda=strrep(labLambda,'.','p');
                ccc = load([foldPD 'S_compPD_lambda_' labLambda improvStr cycleModStr '.mat'],'compPD'); 
                compPD = ccc.compPD;

                %%                
                % prepareToPlot
                oneSt.aSt = prepareToPlot(oneSt.aSt,oneSt.aSt,compPD);
                if plotStackedSimul
                    %% plot
                    % save('./outputs/tmp/preplot_lowerSRT.mat')
                    [joyFigs,oneP.P] = prepareJoyPlot_lowerSRT(oneP.P,oneSt.aSt,aSt_SS_computed,joyFigs);                      
                end
                P.allObj = oneP.P.allObj; P.allObjNames=oneP.P.allObjNames; 
            end
        case 'startupScenario'
                labLambda = num2str(optimal_lambda);

                % load structures containing options
                oneP = load(strcat([P.folderOutput,'P_allLambdas/P_', distrib_at_startupScenario ,'_lambda_',labLambda,improvStr,cycleModStr,'.mat']),'P'); 
                oneP.P = PBEsettings(oneP.P);
                plotsVisible = false;
                [aSt_options,oneP.P] = PBEParams(oneP.P,plotsVisible); 

                % load compPD, used to compute the quartiles
                foldPD = [folderOutput 'S_compPD/']; 
                ccc = load([foldPD 'S_compPD_lambda_' labLambda '.mat'],'compPD'); 
                compPD = ccc.compPD;

                % plot for each lambda
                fprintf('##### Plotting n = %g for lambda %g ##### \n ',i,optimal_lambda)
                % prepareToPlot
                aSt = prepareToPlot(aSt,aSt_options,compPD);

                % load PBE data SS_computed for the same optimal_lambda for joyFigs.fig2 (CFR endpoints)
                folderOutput_SS = strrep(oneP.P.folderOutput,'lowBiom','SS');
                aSt_SS_computed = load(strcat([folderOutput_SS,'aSt_allLambdas/aSt_lambda_',labLambda,'_T_sim',num2str(oneP.P.T_sim),'_nCells',num2str(oneP.P.nCells),cycleModStr,'.mat']),'aSt');
                aSt_SS = aSt_SS_computed.aSt;
                % load compPD, used to compute the quartiles
                foldPD = [folderOutput_SS 'S_compPD/']; 
                ccc = load([foldPD 'S_compPD_lambda_' labLambda '.mat'],'compPD'); 
                compPD_SS = ccc.compPD;
                % prepareToPlot
                aSt_SS = prepareToPlot(aSt_SS,aSt_options,compPD_SS);

                % plot
                if plotStackedSimul
                    joyFigs = prepareJoyPlot(oneP.P,aSt,aSt_SS,joyFigs);

                end
    end
 end
end

if plotSimulation
 % save
 foldFig = './outputs/figures/joyPlot/';
 saveFigs(joyFigs,foldFig)
end


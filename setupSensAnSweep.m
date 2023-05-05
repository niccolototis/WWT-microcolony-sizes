function [P] = setupSensAnSweep(P)
% - Parameter reference
% P.YA = 0.24; % g cell COD formed/ g N oxidized
% P.YH = 0.67; % g cell COD formed/ g COD oxidized
% P.iXB = 0.086; % g N in biomass / g COD in biomass
% P.iXP = 0.06; % g N in endogenous mass/ g COD in endogenous mass
% P.MuH = 6.0; % per day; heterotroph growth
% P.MuA = 0.80; % per day; autotroph growth
% P.bH = 0.62; % per day; heterotroph decay
% P.bA = 0.10; % per day; autotroph decay

% - delete previously opened parallel sessions
% delete(gcp('nocreate'))

P.suppressLegend = false;
P.goodParFound = false; % this can be used in myEvent to stop the ode integration if some trajectory meets some requirements
P.EE.confirmMorrisLikeLHS = false; % I use this to double chech with the profiles what I already see in the surface plot (to doublecheck the cases where with different param combinations i do not see any change)
switch P.parSetup 
 case 'baseline'
    P.nSweeps = 1;  % Number of step changes of the parameters
 case 'chosenSetASM1seq'
    P.nSweeps = 1;  % Number of step changes of the parameters
    P.chosenSetASM1seqName = 'outputs/parameter_modifyOnce/chosenSetASM1seq.xls';
    data = readtable(P.chosenSetASM1seqName);
    goodRows = find(~cellfun(@isnan,num2cell(data.parLims)));
    data = data(goodRows,:);
    P.parToChange = data.Names;
    nVars = length(P.parToChange);       
    P.NewValue = data.NewValue;
    for i = 1:nVars
        str = strcat(['p' num2str(i)]);
        P.(str) = P.NewValue(i);
    end
 case 'paramSweep'
    P.parSweepBaseFileName = 'outputs/parameter_sweep_ASM1seq/ParameterSweep.xls';
    data = readtable(P.parSweepBaseFileName,'Sheet',P.whatSweep); 
    P.pSweep.tableBounds = data;
    P.pSweep.tableBounds.Row = data.Names;
    P.pSweep.parNamesSweep = data.Names;
    nVars = length(P.pSweep.parNamesSweep);
    P.pSweep.k = length(P.pSweep.parNamesSweep);
    P.pSweep.parRef = zeros(1,P.pSweep.k);
    toEval = strcat('P.',data.Names);
    for kk = 1:length(data.Names)
        P.pSweep.parRef(kk) = eval(toEval{kk});            
    end
    P.pSweep.parLims = data.parLims(:)';
    P.pSweep.parLb = data.lb(:)';
    P.pSweep.parUb = data.ub(:)';
    switch P.whatSweep
        case 'LHS'
                P.figure = 'flat';
                P.pSweep.paramsSteps = lhsdesign(P.nSweeps,nVars);
                P.pSweep.paramSetsUsed = nan(size(P.pSweep.paramsSteps));
                for i = 1:nVars
                    str = strcat(['p' num2str(i)]);
                    P.pSweep.paramSetsUsed(:,i) = P.pSweep.parLb(i)+P.pSweep.paramsSteps(:,i)*(P.pSweep.parUb(i)-P.pSweep.parLb(i));
%                     P.pSweep.paramsSteps(:,i) = P.pSweep.parRef(i)+P.pSweep.paramsSteps(:,i)*(P.pSweep.parLims(i)-P.pSweep.parRef(i));
                    P.(str) = P.pSweep.paramSetsUsed(:,i);
                end
        case 'Ladder'
                P.nSweeps = 20;
                oneColSteps = linspace(0,1,P.nSweeps+1)';
                oneColSteps = oneColSteps(2:end);
                P.pSweep.paramsSteps = repmat(oneColSteps,1,nVars);
                for i = 1:nVars
                    str = strcat(['p' num2str(i)]);
                    P.pSweep.paramsSteps(:,i) = P.pSweep.parRef(i)+P.pSweep.paramsSteps(:,i)*(P.pSweep.parLims(i)-P.pSweep.parRef(i));
                    P.(str) = P.pSweep.paramsSteps(:,i);
                end
        case 'meshgrid'
            % Set up Parameter Sweep with meshgrid
            % - Create a 2-D grid of parameters by using the meshgrid function.
            p1 = sort(linspace(parRef(1), parLims(1), P.nSweeps));
            p2 = sort(linspace(parRef(2), parLims(2), P.nSweeps));
            [p1,p2] = meshgrid(p1,p2);
            P.p1 = p1;
            P.p2 = p2;
    end   
    P.nRuns = numel(P.p1); 
    
    % - General operations for paramSweep option
    if length(P.pSweep.parNamesSweep)~= 2 && strcmp(P.figure,'surface')
        ME = MException('Need to have 2 parameter to plot the surface plot');
        throw(ME);
    end
    
    % define plot title
    P.title = [];
    for jj = 1:length(P.pSweep.parNamesSweep) % for each of the parameters that I change
        E = strcat(P.pSweep.parNamesSweep{jj},'=');
        E = strcat([' ',E,sprintf('%.4g', P.pSweep.parRef(jj))]); %I assign the ii-th par value to its specific parameter name 
        P.title = [P.title E];
    end
    disp(P.title)
    
    % Still load the data of Morris for the lb and ub on biomass
    % concentration (used in runASM1seqparamSweep.m)
    
    P.outputMorrisFileName = './outputs/parameter_Morris/Parameters_Morris.xls';
    data = readtable(P.outputMorrisFileName,'Sheet','runs','Range','A:D');        
    goodRows = find(~cellfun(@isempty,data.Names));
    data = data(goodRows,:);
    data.Row = data.Names;
    P.pSweep.lb_biom_conc = data.lb('biom_conc');
    P.pSweep.ub_biom_conc = data.ub('biom_conc');
    P.pSweep.ref_biom_conc = data.RefValue('biom_conc');
    
 case 'EEMorris'
    P.figure = 'surface';     
    % INPUTS:
    % - lb: lower bound on parameter values
    % - ub: upper bound on parameter values
    % - M: number of Morris trajectories generated
    % - p: number of grid points
    % - r: number of chosen Morris trajectories

    P.inputMorrisFileName = './inputs/parameter_Morris/Parameters_Morris.xls';
    P.outputMorrisFileName = './outputs/parameter_Morris/Parameters_Morris.xls';
    if ~exist(P.outputMorrisFileName,'file')
        copyfile(P.inputMorrisFileName,P.outputMorrisFileName)
    end
    
    data = readtable(P.outputMorrisFileName,'Sheet','runs');
    
    % - Clear previously run simulations
    clearFrom = find(strcmp(data.Properties.VariableNames,'nRun'));
    clearEnd = length(data.Properties.VariableNames);
    clarCols = clearEnd-clearFrom;
    tableToClear = array2table(strings(size(data,1)+1,clarCols));
    writetable(tableToClear,P.outputMorrisFileName,'WriteRowNames',false,'WriteVariableNames',false,'Sheet','runs','Range','G1');
    
    % - Take relevant data
    goodRows = find(~cellfun(@isempty,data.Names));
    data = data(goodRows,:);
    P.EE.parNamesMorris = data.Names;
    P.EE.k = length(P.EE.parNamesMorris);
    P.parRef = zeros(1,P.EE.k);
    toEval = strcat('P.',data.Names);
    for kk = 1:length(data.Names)
        P.parRef(kk) = eval(toEval{kk});            
    end
    % ub and lb need to be given as row vectors
    P.EE.lb = data.lb(:)';
    P.EE.ub = data.ub(:)';
    P.EE.p = 10; % number of levels each parameter
    
    % - ﻿We start by generating a high number of different Morris trajectories, M~500-1000
    P.EE.M = 30; % All the trajectories generated 
    
    % -﻿The EE method is based on the construction of r trajectories in the input space, typically between 10 and 50
    P.EE.r = 20; % number of traj selected   
    
    P.nRuns = P.EE.r*(P.EE.k+1);            
    
    filename_comb = strcat(['outputs/parameter_Morris/combMorris_k',num2str(P.EE.k),'_p',num2str(P.EE.p),'_M',num2str(P.EE.M),'_r',num2str(P.EE.r),'.mat']);
    if isfile(filename_comb)
        load(filename_comb);
    else
        comb = nchoosek(1:P.EE.M,P.EE.r);
        save(filename_comb,'comb','-v7.3');
    end      
    
    P = runMorris(P,comb);
end

P = defineColorMap(P);
P.goodParFound = false;
P.thisLev = 1; % Counter to iterate inside while loop in runASM1seqModel_ASM1
end




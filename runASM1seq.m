function [P,varargout] = runASM1seq(P,varargin)
global gcs folderOutput T_sim cycleModStr improvStr 
% Initialize options, counters and flags
P = addToRecover(P);
P = settingsASM1seq(P); 

if P.plotASM1seq && ~P.plotPhases
 figHET = figure('visible','on','Position',[892,385,560,420]);
 figHET.Name = 'figHET';
 figAUT = figure('visible','on','Position',[881,54,560,420]);
 figAUT.Name = 'figAUT';
end

if P.runASM1seq 
 % If [true,false] this fo  r loop runs the whole bag model twice to make sure that the
 % two implementations with state vector representing A)weigh and
 % B)concentrations produce the same trajectories
 P.gc = char(gcs(P.variableInGrams));   % string to create and access the structure fields    
 P = defineKineticParams(P);
 % - In setupSensAnSweep_parallel some of the kineticParams are changed
 % plus lead to changes in initializeSmatrix and defineInletConcentrations,
 % which are inside runASM1seqForOneParam  
 P = setupSensAnSweep(P); 
 
 switch P.parSetup
    case 'baseline'
        [t_res_base,X_c_res_base,S_c_res_base,V_res_base] = runASM1seqbaselineParams(P,'baseline');
        % - organize output
        varargout{1} = t_res_base;varargout{2}=X_c_res_base;varargout{3}=S_c_res_base;varargout{4}=V_res_base;
        
    case 'chosenSetASM1seq'
        P = changeParam(P,1,P.parToChange);
        [t_res,X_c_res,S_c_res,V_res] = runASM1seqForOneParam(1,P);  
        
        % - write to file   
        fName = [folderOutput,'resASM1seq_',num2str(T_sim),improvStr,cycleModStr,'.mat'];            
        fName = checkNotOverWrite(fName);
        
        save(fName,'t_res','X_c_res','S_c_res','V_res', '-v7.3');

        % - organize output
        varargout{1} = t_res;varargout{2}=X_c_res;varargout{3}=S_c_res;varargout{4}=V_res;

        % - add line with new parameters to plot
        if P.plotASM1seq
            P = plotOneRunASM1seq(P,figHET,X_c_res(:,P.HET),'HET','Chosen set',P.parSetup,[]);
            P = plotOneRunASM1seq(P,figAUT,X_c_res(:,P.AUT),'AUT','Chosen set',P.parSetup,[]);
            adjustPlotAndSave(P,figHET,'HET');
            adjustPlotAndSave(P,figAUT,'AUT');                
        end
 case 'paramSweep'
    runASM1seqparamSweep(P);
    varargout{1} = [];varargout{2}=[];varargout{3}=[];varargout{4}=[];
    
 case 'EEMorris'      
    %         % - try just one run to debug
    %         for runNo = 1:42
    %             oneResult = runASM1seqForOneParam(runNo,P);
    %         end
 
    [P,modelOutputs] = parallelRuns_ASM1seq(P);        

    % - Compute Elementary Effects 
    % - In EEtot_thisPRange I consider
    % that changes in parameters within each range as having the same
    % kind of importance, so if the effect on the output is the same
    % for 2 parameters, even if they have different ranges (1 narrow, 1
    % large) I want their effect to be the same here
    [EEtot] = deal(zeros(P.EE.k,P.EE.nTraj));     
    for z = 1:P.EE.nTraj
        EE = zeros(1,P.EE.k); % row vector
        for i = 2:(P.EE.k+1) % along the rows
            ind = find(P.EE.X(i,:,z)~=P.EE.X(i-1,:,z)); %@NT I am finding the column index that differentiates this row from the previous one
            EE(ind) = (modelOutputs(i,z)-modelOutputs(i-1,z))/(P.EE.X_Params(i,ind,z)-P.EE.X_Params(i-1,ind,z)); %@NT ind will take all values between 1 and k in random order
        end
        EEtot(:,z) = EE;
    end
    parRange = P.EE.ub-P.EE.lb;
    EEtot_thisPRange = EEtot.*parRange(:);
    
    % Calculate mu and sigma for each input parameter
    % Note that the absolute value of the elementary effects is taken while
    % calculating mu!
    [failedEEs,mu,mu_abs,sigma,mu_thisRange,mu_abs_thisRange,sigma_thisRange] = deal(zeros(size(EEtot,1),1));
    for i = 1:size(EEtot,1)
        oneRow = EEtot(i,:); 
        % remove NaN, take into consideration only good runs for the computation of EEs
        failedEEs(i) = sum(isnan(oneRow));
        oneRow = oneRow(~isnan(oneRow));
        mu(i) = mean(oneRow);
        mu_abs(i) = mean(abs(oneRow));
        sigma(i) = std(oneRow);
        
        % - EEtot_thisPRange NB!
        oneRow = EEtot_thisPRange(i,:);
        % remove NaN, take into consideration only good runs for the computation of EEs
        oneRow = oneRow(~isnan(oneRow));
        mu_thisRange(i) = mean(oneRow);
        mu_abs_thisRange(i) = mean(abs(oneRow));
        sigma_thisRange(i) = std(oneRow);
    end
    outputTable = array2table([failedEEs,mu,mu_abs,sigma,mu_thisRange,mu_abs_thisRange,sigma_thisRange],'VariableNames',{'failedEEs','mu','mu_abs','sigma','mu_thisRange','mu_abs_thisRange','sigma_thisRange'});
    save('outputs/parameter_Morris/EE_computed.mat','outputTable');
    writetable(outputTable,P.outputMorrisFileName,'Sheet','EE','Range','F1')        
 end % end switch
else % If ~runASM1seq
 % - Run the ancillary scripts to have global variables defined for later usage
 P = settingsASM1seq(P); 
 P = defineKineticParams(P);
 % - I run setupSensAnSweep one time before initial
 % conditions because I need to have P.paramSetup defined. In
 % alternative put the definition of P.parSetup in generalSettings
 if strcmp(P.parSetup,'bestSweep_ASM1seq_to_PBE')
    [P,bestSweep_ASM1seq] = setupSensAnSweep(P); 
    P = initializeSmatrix(P);
    P = defineInletConcentrations(P);
    P = defineInitialConditions_t0(P,bestSweep_ASM1seq);    
 else       
    P = setupSensAnSweep(P);
    P = initializeSmatrix(P); 
    P = defineInletConcentrations(P);
    P = defineInitialConditions_t0(P);
 end
 
 % - Load previously saved data
 % - Assign parameters for this run
 switch P.parSetup
    case 'bestSweep_ASM1seq_to_PBE'
        vv = ['v' num2str(P.nChosenParSet_ASM1seq) '_'];
        load([folderOutput,'resASM1seq_alternative_',vv,num2str(T_sim),improvStr,cycleModStr,'.mat'],'t_res','X_c_res','S_c_res','V_res');           
    case 'chosenSetASM1seq'
        P = changeParam(P,1,P.parToChange);
        if ~isfile([folderOutput,'resASM1seq_',num2str(T_sim),improvStr,cycleModStr,'.mat'])
            error('Run the ASM1seq model first')
        else
            load([folderOutput,'resASM1seq_',num2str(T_sim),improvStr,cycleModStr,'.mat'],'t_res','X_c_res','S_c_res','V_res');
        end
        if P.plotASM1seq
            if P.plotPhases
                plotOnlyPhases(P);
                plotSSforPaper(P,X_c_res,{'HET' 'AUT'});
            else
                P = plotOneRunASM1seq(P,figHET,X_c_res(:,P.HET),'HET','Chosen set','Chosen set',1);
                P = plotOneRunASM1seq(P,figAUT,X_c_res(:,P.AUT),'AUT','Chosen set','Chosen set',1);
                adjustPlotAndSave(P,figHET,'HET');
                adjustPlotAndSave(P,figAUT,'AUT');
            end
        end
    case 'baseline'
        load([folderOutput,'resASM1seq_BASE_',num2str(T_sim),improvStr,cycleModStr,'.mat'],'t_res_base','X_c_res_base','S_c_res_base','V_res_base');
        t_res = t_res_base;
        X_c_res = X_c_res_base;
        S_c_res = S_c_res_base;
        V_res = V_res_base;
        if P.plotASM1seq
            P = plotOneRunASM1seq(P,figHET,X_c_res(:,P.HET),'HET','Nominal - literature','Nominal - literature',1);
            P = plotOneRunASM1seq(P,figAUT,X_c_res(:,P.AUT),'AUT','Nominal - literature','Nominal - literature',1);
            adjustPlotAndSave(P,figHET,'HET');
            adjustPlotAndSave(P,figAUT,'AUT');                
        end
    case 'paramSweep'
        runASM1seqparamSweep(P);
        varargout{1} = [];varargout{2}=[];varargout{3}=[];varargout{4}=[];
    otherwise
        error('Select an appropriate option for parSetup')
 end

 % - assign output
 varargout{1} = t_res;varargout{2}=X_c_res;varargout{3}=S_c_res;varargout{4}=V_res;
end
return

function P = addToRecover(P)
% Here I assign the fields of P to make this file able to call the
% dependent funtions present in run..._bis
P.modelStructure = 'actualSystem'; 
P.modelICPars = 'actualSystem';
P.modelKinPars = 'actualSystem';
P.ICisSS = true; 
P.ICisSCells = false;
return

function fName = checkNotOverWrite(fName)
global improvStr cycleModStr iT parametersEffectOnMSD
 % If I am running a new batch of experiments the same day, I
 % want to avoid overwriting the results = > I shift the index iT
 % to the next one available
 if parametersEffectOnMSD 
    while exist(fName,'file')
        iT = iT+1;
        improvStr = ['_impro_' char(datetime('today')) '_v' num2str(iT)];
        fName = [folderOutput,'resASM1seq_',num2str(T_sim),improvStr,cycleModStr,'.mat'];
    end
 end
return








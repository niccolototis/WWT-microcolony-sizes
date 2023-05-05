function varargout = runASM1seqForOneParam(aRun,P,varargin)
global myPlotOff nnsos step 
start_time = tic;
fprintf('Start integration of run %g \n',aRun)

if nargin ==5
 X_c_tall = varargin{1};
 S_c_tall = varargin{2};
 V_tall = varargin{3};
end

% - Assign parameters for this run
switch P.parSetup
 case 'chosenSetASM1seq'
    % - do nothing this should already be done in
    % runASM1seq.m . In this case I do it in runASM1seq to
    % have that parameters are updated inside struct P -> the correct values go to PBE 
 case 'paramSweep'
    for jj = 1:length(P.pSweep.parNamesSweep) % for each of the parameters that I change
        str = strcat(['p' num2str(jj)]);       
        P.(char(P.pSweep.parNamesSweep(jj))) = P.(str)(aRun); %I assign the ii-th par value to its specific parameter name 
    end
 case 'EEMorris'
    % - following indexing is correct
    [ind_x,ind_z] = ind2sub([P.EE.k+1 P.EE.nTraj],aRun);
    theseParams = P.EE.X_Params(ind_x,:,ind_z);
    for jj = 1:length(P.EE.parNamesMorris) % for each of the parameters that I change        
        P.(char(P.EE.parNamesMorris(jj))) = theseParams(jj); %I assign the ii-th par value to its specific parameter name 
    end
end

% - The following functions make use of the parameters assigned for this
% specific run
P = initializeSmatrix(P); 
P = defineInletConcentrations(P);
% Define initial conditions (at time zero of the first phase)
P = defineInitialConditions_t0(P);


% - Initialize at the start of the first cycle
t_all = []; N_wbag_tall=[]; 
P.indPhase = 1; 
options = odeset('NonNegative',true,'RelTol', 1e-10, 'AbsTol', 1e-10, 'MaxSTep', 1000);

while P.indPhase<= P.nPhases_all 

 % - Update control for this phase
 P = controlOperations(P);
 
 % - Prepare initial conditions for the upcoming phase
 P = defineInitialConditions_intermediateSteps(P);
 P = checksInterm(P);
 P.tspanThisFrag = P.tspans_ASM1seq.(strcat(['f',num2str(P.indPhase)]));
 
 % - plot output to console only when I am not in parallel mode
 % otherwise it gets messy
 if P.runningParallel ==false
    disp(strcat(['time = ',num2str(P.tspanThisFrag(1)),' phase=',num2str(P.phaseCtrl), ' onAer=',num2str(P.onAer),' o2conc',num2str(P.O2_c_ref),' pastPhase:',num2str(P.pastPhase)]))
 end

 step = 0;
 % - Run ODEs
 [t_oneFrag, N_wbag_oneFrag] = ode15s(@(t, x_all) odeFunASM1seq(t, x_all,P),P.tspanThisFrag,P.this_x_all_0,options); 

 % - Store results
 t_all = [t_all;t_oneFrag];
 N_wbag_tall = [N_wbag_tall;N_wbag_oneFrag];
 P.this_x_all_0 = N_wbag_oneFrag(end,:)';
 P.indPhase = P.indPhase+1;
end    

% -- Organize output
if (P.nPart+P.nSolub+1)~= size(N_wbag_tall,2)
 ME = MException('Check nPart and nSolub indexes');
 throw(ME);
end
V_tall_gc = N_wbag_tall(:,end); % The last column is volume keeps unchanged
X_tall = N_wbag_tall(:,1:P.nPart);
S_tall = N_wbag_tall(:,(P.nPart+1):(P.nPart+P.nSolub));

if P.variableInGrams 
 % convert to concentrations
 X_c_tall_gc = X_tall./V_tall_gc;                                                   
 S_c_tall_gc = S_tall./V_tall_gc;
else
 % do nothing
 X_c_tall_gc = X_tall;
 S_c_tall_gc = S_tall;
end


switch P.parSetup
 case {'chosenSetASM1seq','baseline'}
    % -- Print to console 
    rqrmts_gc = zeros(length(t_all),length(P.rqrmNames));
    rqrmts_gc(:,P.tCOD) = S_c_tall_gc(:,P.sI) + S_c_tall_gc(:,P.sS) ;
    rqrmts_gc(:,P.tN) = S_c_tall_gc(:,P.sNH4) + S_c_tall_gc(:,P.sNO) + S_c_tall_gc(:,P.sND);
    rqrmts_gc(:,P.MLVSS) = (X_c_tall_gc(:,P.xI) + X_c_tall_gc(:,P.xS) + X_c_tall_gc(:,P.HET) + X_c_tall_gc(:,P.AUT) + X_c_tall_gc(:,P.xP))/1.42; 
    rqrmts_gc(:,P.OUR) = -1*S_c_tall_gc(:,P.sO2);
    disp(strcat(['at the steady state the tCOD = ',num2str(rqrmts_gc(end,P.tCOD)),' mg/L  ,  totalN=',num2str(rqrmts_gc(end,P.tN)),' mg/L']))
    disp(strcat(['Requirements tCOD<',num2str(125),' mg/L  ,  totalN<',num2str(15),' mg/L']))
    % -- Plot
    if ~myPlotOff
        X_c_tall.g = X_c_tall_gc; S_c_tall.g=S_c_tall_gc; V_tall.g=V_tall_gc;
        nnsos = 1;
        plotASM1seqresults(P,t_all,X_c_tall,S_c_tall,V_tall)
        nnsos = [1 2];
    end
    
     % -- Organize output for the single run --> to PBE
    varargout{1} = t_all;varargout{2}=X_c_tall_gc;varargout{3}=S_c_tall_gc;varargout{4}=V_tall_gc;
    
 case 'EEMorris'
    diffAUT_fromSS = abs((X_c_tall_gc(end,P.AUT)-X_c_tall_gc(1,P.AUT))/X_c_tall_gc(1,P.AUT));
    diffHET_fromSS = abs((X_c_tall_gc(end,P.HET)-X_c_tall_gc(1,P.HET))/X_c_tall_gc(1,P.HET));
    diffTOT = diffAUT_fromSS+diffHET_fromSS;
    diffTOT = -diffTOT; % I make it negative because I want to have the positive effect of a parameter to be the reduction of diffAUT_fromSS+diffHET_fromSS
    oneResult = diffTOT;
    % -- Organize output 
    varargout{1} = oneResult; varargout{2}=[]; varargout{3}=[]; varargout{4}=[]; varargout{5}=[]; 

 case 'paramSweep'
    switch P.figure
        case 'surface'
            oneResult = (X_c_tall_gc(end,P.AUT)-X_c_tall_gc(1,P.AUT))*(100/X_c_tall_gc(1,P.AUT));
        case 'flat'
            switch P.plot
                case 'AUT'
                    oneResult = X_c_tall_gc(:,P.AUT);
                case 'HET'
                    oneResult = X_c_tall_gc(:,P.HET);
                case 'both_HET_AUT'
                    oneResult = [X_c_tall_gc(:,P.HET);X_c_tall_gc(:,P.AUT)];
                case 'nitrogen'
            end
    end
    
    % -- Organize output 
    X_c_end = X_c_tall_gc(end,:);
    S_c_end = S_c_tall_gc(end,:);

    varargout{1} = oneResult; varargout{2}=X_c_tall_gc; varargout{3}=S_c_tall_gc; varargout{4}=X_c_end; varargout{5}=S_c_end; 
end

% Print end statement
end_time = toc(start_time);
fprintf('Run %g completed in %g seconds. \n',aRun,end_time)
end

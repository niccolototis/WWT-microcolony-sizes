function aSt = runPBEforOneParam(aSt,P,M_N_TI_0_full,t_all_ASM1seq,mu_ASM1seq_tall,Ker_B,Ker_V,sp) 
global step steP IC_at distrib_at_startupScenario cycleModStr improvStr
 
% - Initialize data storage matrices
N_TI_res_from_ICnorm = zeros(P.ntpts_ASM1seq,P.nCells*P.nChar); % 
[x1_all_res, w1_all_res] = deal(zeros(P.ntpts_ASM1seq,P.nChar));
aSt.(sp).x2_all_res = repmat(aSt.(sp).x2_all',P.ntpts_ASM1seq,1);

for aChar = 1:P.nChar
 % 4Loop . take the N~ for each characteristic curve for all species. NB taking the fist
 % characteristics does not mean that x1 and w1 are the same for all curves!
 N_TI_x_w_t0_aChar = M_N_TI_0_full.(sp)(:,aChar); % slice the matrix taking just one characteristics
 options = odeset('NonNegative',1,'RelTol',1e-6,'AbsTol', 1e-6);
 t_all_PBE = zeros(P.ntpts_ASM1seq,1); N_res_aChar=zeros(P.ntpts_ASM1seq,length(N_TI_x_w_t0_aChar)); 
 %---------
 P.indPhase = 1; P.thisTpoint=1; P.broken=false; steP=0;
 while P.indPhase<= P.nPhases_all         
    step = 0;
    
    % -- update control for this phase
    tspanThisFrag = P.tspans_ASM1seq.(strcat(['f',num2str(P.indPhase)]));
    P = controlOperations(P);
    disp(strcat(['time = ',num2str(tspanThisFrag(1)),' phase=',num2str(P.phaseCtrl), ', pastPhase: ',P.pastPhase,' - thisPhase: ',P.thisPhase]))
    
    % -- prepare initial conditions for this phase
    [N_TI_x_w_t0_aChar,P] = defineInitialConditions_intermediateSteps_PBE(N_TI_x_w_t0_aChar,P);
    
    %---------
    % Change this particular warning into an error
    % warnId = 'MATLAB:ode15s:IntegrationTolNotMet';
    % warnstate = warning('error', warnId);
    % That can be caught with try/catch
    % try  % Quick-failing integration, suggested by horchler
    [t_oneFrag, N_res_aChar_oneFrag] = ode15s(@(t,N_all_aChar) odeFunPBE(t,N_all_aChar,aSt.(sp),Ker_B.(sp),Ker_V.(sp),t_all_ASM1seq,mu_ASM1seq_tall,aChar,sp,P), tspanThisFrag, N_TI_x_w_t0_aChar,options); 
    %         catch ME
    %             % Check if we indeed failed on meeting the tolerances
    %             disp(strcat(['Error during integration at step = ',num2str(steP)]))
    %             if strcmp(ME.identifier, warnId)
    %                 disp('Remember that Subs_V_whole_tall are discontinuous functions that are fed to the ode system. You can try reduce dt_ASM1seq to have it smoother')
    %             else
    %                 throw(ME);
    %             end
    %         end
    %         warning(warnstate); % Don't forget to reset the warning status

    % -- store results
    positions = (P.tFragm_lengths(P.indPhase)+1):P.tFragm_lengths(P.indPhase+1);
    t_all_PBE(positions) = t_oneFrag;
    N_res_aChar(positions,:) = N_res_aChar_oneFrag;
    N_TI_x_w_t0_aChar = N_res_aChar_oneFrag(end,:)';
    N_TI_x_w_t0_aChar = fixLowerThreshold_fast(N_TI_x_w_t0_aChar,P);
    P.indPhase = P.indPhase+1;
 end  

 % -- Organize output
 % N_TI_res_from_ICnorm the name is explanatory: the from_ICnorm means that the PBE is initialized 
 % with an initial distribution which is a pdf. However we are not
 % using the normalized PBE, thus the output has to be considered as a
 % NDF (If there is particle proliferation the zeroth moment grows above 1)
 N_TI_res_from_ICnorm(:,((aChar-1)*(P.nCells)+(1:P.nCells))) = N_res_aChar(:,(1:P.nCells));
 x1_all_res(:,aChar) = N_res_aChar(:,(P.nCells+1)); 
 w1_all_res(:,aChar) = N_res_aChar(:,end);
end

% -
if P.tspanCut    
 % - adjust P.idx_tpts_plot_PBE for when I compute moments
 P.tmp_idx_tpts_plot_PBE = [1 P.idx_tpts_plot_PBE(t_all_PBE(P.idx_tpts_plot_PBE)~=0)];    
 % - recreate the full time vector, which otherwise would have zeros
 % after the specified tspanCut until t = 50 (ss)   
 indPhase_tmp = 1;
 while indPhase_tmp<= P.nPhases_all_ASM1seq         
    tspanThisFrag = P.tspans_ASM1seq.(strcat(['f',num2str(indPhase_tmp)]));
    positions = (P.tFragm_lengths(indPhase_tmp)+1):P.tFragm_lengths(indPhase_tmp+1);
    t_all_PBE(positions) = tspanThisFrag;
    indPhase_tmp = indPhase_tmp+1;
 end
end
% -
[t_all_PBE,N_TI_res_from_ICnorm,x1_all_res,w1_all_res] = shrinkSize(P,t_all_PBE,N_TI_res_from_ICnorm,x1_all_res,w1_all_res);
if ~isequal(string(N_TI_res_from_ICnorm(1,:)'),string(P.M_N_TI_0))
 ME = MException('should be the same, check if shrinkSize is removing the first time point (initial conditions)');
 throw(ME);
end
aSt.(sp).t_all_PBE = t_all_PBE;
aSt.(sp).N_TI_res_from_ICnorm = N_TI_res_from_ICnorm;
aSt.(sp).x1_all_res = x1_all_res;
aSt.(sp).w1_all_res = w1_all_res;
% -
[m1_x2,m2_x2,Marg_N_x1_RESC0to1,Marg_N_x2_RESC0to1,N_tot_t_all,Biomass_tot] = storePBEOutput(P,N_TI_res_from_ICnorm,x1_all_res,w1_all_res,aSt,sp);
aSt.(sp).m1_x2 = m1_x2;
aSt.(sp).m2_x2 = m2_x2;
aSt.(sp).Marg_N_x1_RESC0to1 = Marg_N_x1_RESC0to1;
aSt.(sp).Marg_N_x2_RESC0to1 = Marg_N_x2_RESC0to1;
aSt.(sp).Biomass_tot = Biomass_tot;
aSt.(sp).N_tot_t_all = N_tot_t_all;
% -
aSt = from_N_to_pdf_pvf(aSt,sp);

if strcmp(sp,P.spNames(P.speciesPBE(end))) 
% if P.saveOptimalRun
%     save(strcat([P.folderData,'aSt_',P.thisObjFunType,'_lambda',num2str(P.lambda),'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells),cycleModStr,'.mat']),'aSt');
% else
 if ~exist([P.folderOutput,'P_allLambdas/'],'dir')
    mkdir([P.folderOutput,'P_allLambdas/'])
 end
 if ~exist([P.folderOutput,'aSt_allLambdas/'],'dir')
    mkdir([P.folderOutput,'aSt_allLambdas/'])
 end
 if(strcmp(IC_at,'startupScenario'))
    save(strcat([P.folderOutput,'P_allLambdas/P_', distrib_at_startupScenario ,'_lambda',num2str(P.lambda),improvStr,cycleModStr,'.mat']),'P'); 
    save(strcat([P.folderOutput,'aSt_allLambdas/aSt_', distrib_at_startupScenario ,'_lambda',num2str(P.lambda),'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells),improvStr,cycleModStr,'.mat']),'aSt');
 else
    save(strcat([P.folderOutput,'P_allLambdas/P_lambda_',num2str(P.lambda),improvStr,cycleModStr,'.mat']),'P'); 
    save(strcat([P.folderOutput,'aSt_allLambdas/aSt_lambda_',num2str(P.lambda),'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells),improvStr,cycleModStr,'.mat']),'aSt');
 end
end
return

function [N_TI_x_w_t0_aChar,P] = defineInitialConditions_intermediateSteps_PBE(N_TI_x_w_t0_aChar,P)
N_TI_aChar = N_TI_x_w_t0_aChar(1:P.nCells);
if P.instantSludgeRemoval % I washout particles instantanously
 if (strcmp(P.modelWaterSludgeSegreg,'segregated') && P.onPurge) || (strcmp(P.modelWaterSludgeSegreg,'stirred') && P.onEfflux) % @NB in this way the purge is implemented at the beginning of the purge window 
    Dw = N_TI_aChar*(1/P.SRT)*(1/P.cyclesOneDay); 
    N_TI_aChar = N_TI_aChar-Dw; % Death does not affect the redistribution among different cells in the CAT
 end  
else
 P.Dw_purge = 0;
 if strcmp(P.modelWaterSludgeSegreg,'segregated') && P.onPurge
    %- remove 1/SRT for each day, so that on average I have full
    %sludge renewal in a SRT days period
    P.Dw_purge = N_TI_aChar*(1/P.SRT)/P.purge_tDays;   
 else
    % do nothing, the "stirred option is implemented directly inside the ode function"
 end
end
N_TI_x_w_t0_aChar(1:P.nCells) = N_TI_aChar;
return

function [t_all_PBE,N_TI_res_from_ICnorm,x1_all_res,w1_all_res] = shrinkSize(P,t_all_PBE,N_TI_res_from_ICnorm,x1_all_res,w1_all_res)
 t_all_PBE = t_all_PBE(P.idx_tpts_plot_PBE);
 N_TI_res_from_ICnorm = N_TI_res_from_ICnorm(P.idx_tpts_plot_PBE,:);
 x1_all_res = x1_all_res(P.idx_tpts_plot_PBE,:);
 w1_all_res = w1_all_res(P.idx_tpts_plot_PBE,:);
return

function [m1_x2_RESC0to1,m2_x2_RESC0to1,Marg_N_x1_RESC0to1,Marg_N_x2_RESC0to1,N_tot_t_all,Biomass_tot] = storePBEOutput(P,N_TI_res_from_ICnorm,x1_all_res,w1_all_res,aSt,sp)
% -- Store output

% - initialize vectors and matrices
 [M_N_res_ICnorm,M_N_TI_res_ICnorm,M_N_res_RESC0to1_ICnorm,M_N_TI_res_RESC0to1_ICnorm] = deal(zeros(P.nCells,P.nChar,P.ntpts_plot_PBE));
 Marg_N_x1_RESC0to1 = zeros(P.ntpts_plot_PBE,P.nChar);
 Marg_N_x2_RESC0to1 = zeros(P.ntpts_plot_PBE,P.nCells);
 [m1_x2_RESC0to1,m2_x2_RESC0to1,N_tot_ICnorm,m0_RESC0to1] = deal(zeros(P.ntpts_plot_PBE,1));

% Matrix organization of N = M_N and Ñ=M_N_TI 
%  __________
% |\         \
% | \ntpts_plot_PBE        
% |  \         \
% |   \ __nChar_\
%  \   |        |
%   \  |nCells  |
%    \ |        |
%     \|________|

if P.tspanCut
 P.ntpts_plot_PBE = length(P.tmp_idx_tpts_plot_PBE);
end

for jj = 1:P.ntpts_plot_PBE
 P.M_NisTilde = true;
 % - Reshape into matrix
 M_N_TI_res_ICnorm(:,:,jj) = reshape(N_TI_res_from_ICnorm(jj,:),P.nCells,P.nChar); % for every time point I reconstruct a matrix for N tilde

 % - Convert Ñ -> N 
 M_N_res_ICnorm(:,:,jj) = M_N_TI_res_ICnorm(:,:,jj)./w1_all_res(jj,:); % converting the matrix Ntilde to the matrix of N       

 % - Compute zeroth_moment = Ntot for each time point. This N_tot_from_ICnorm is the Ntot when initial condition a pdf => at t=0: zeroth_moment=Ntot=1
 N_tot_ICnorm(jj) = calcMoment(M_N_TI_res_ICnorm(:,:,jj),P.nChar);

 % - Get the normalized N_norm and Ñ_norm (memo this is obtained with a stardard PBE starting from normalized IC (at t = 0 ndf is a pdf with ntot_0=1) and then dividing ndf by ntot(t)) 
 M_N_TI_res_RESC0to1_ICnorm(:,:,jj) = M_N_TI_res_ICnorm(:,:,jj)/N_tot_ICnorm(jj);
 M_N_res_RESC0to1_ICnorm(:,:,jj) = M_N_res_ICnorm(:,:,jj)/N_tot_ICnorm(jj);

 % Compute first moment ( = E[x1]) and second moment of x1
 % m1_x1_TI_RESC0to1_allChar = sum(M_N_TI_res_normalized_ICnorm(:,:,jj).*aSt.(sp).x1_all,1);
 % m2_x1_TI_RESC0to1_allChar = sum(M_N_TI_res_normalized_ICnorm(:,:,jj).*aSt.(sp).x1_all.^2,1);
 % % - Integral across the columns (x1 coord)
 % M_NisTilde = true;
 % m1_x1_RESC0to1(jj) = computeMarg(m1_x1_TI_RESC0to1_allChar,M_NisTilde,[],[],nChar,'x2');
 % m2_x1_RESC0to1(jj) = computeMarg(m2_x1_TI_RESC0to1_allChar,M_NisTilde,[],[],nChar,'x2');

 % - Compute first moment ( = E[x2]) and second moment of x2
 m0_RESC0to1(jj) = calcMoment(M_N_TI_res_RESC0to1_ICnorm(:,:,jj),P.nChar);
 reallyEqual(m0_RESC0to1(jj),1,'~','The zeroth moment of the normalized ndf should be 1 at each time point');
 m1_x2_RESC0to1(jj) = calcMoment(M_N_TI_res_RESC0to1_ICnorm(:,:,jj).*aSt.(sp).x2_all,P.nChar);
 m2_x2_RESC0to1(jj) = calcMoment(M_N_TI_res_RESC0to1_ICnorm(:,:,jj).*aSt.(sp).x2_all.^2,P.nChar);

 % - Compute marginals of N_norm
 Marg_N_x1_RESC0to1(jj,:) = sum(M_N_res_RESC0to1_ICnorm(:,:,jj),1); % - Here I consider N 
 Marg_N_x2_RESC0to1(jj,:) = computeMarg(M_N_TI_res_RESC0to1_ICnorm(:,:,jj),P.M_NisTilde,[],[],P.nChar,'x2'); % - Here I consider Ñ 
    
 % - Compute the coefficient of variation ?
 % CV2_x2 = sqrt(m2_x2-(m1_x2).^2)./m1_x2;
end

% - Checks
if P.nChar ==1
 jj = 1;
 reallyEqual(M_N_res_ICnorm(:,:,jj),P.M_N_0,' = ','should be the same, check if shrinkSize is removing the first time point (initial conditions)');
 reallyEqual(M_N_TI_res_ICnorm(:,:,jj),P.M_N_TI_0,' = ','should be the same, check if shrinkSize is removing the first time point (initial conditions)');
 reallyEqual(M_N_res_RESC0to1_ICnorm(:,:,jj),P.M_N_0,' = ','should be the same, check if shrinkSize is removing the first time point (initial conditions)');
 reallyEqual(M_N_TI_res_RESC0to1_ICnorm(:,:,jj),P.M_N_TI_0,' = ','should be the same, check if shrinkSize is removing the first time point (initial conditions)');
 reallyEqual(Marg_N_x2_RESC0to1(jj,:)',P.M_N_0,' = ','should be the same, check if shrinkSize is removing the first time point (initial conditions)');
end
intMarg = sum(Marg_N_x2_RESC0to1,2);
diffs = intMarg-ones(size(intMarg));
if any(diffs>1e-13)
 ME = MException('Marg_N_x2_RESC0to1 (zeroth moments in all the cells in the grid) sums to 1, like the pdf integrates to 1');
 throw(ME);
end
% -- Print inconsistencies found 
printInconsistencies(P); % @MEMO need to fix this so that it prints properly

% - Compute the total change in size. Note, importantly I have 2 normalization constants: 
% 1 - N_tot_from_ICnorm (Increment of Ntot through the integration of the PBE, started from IC: ndf = pdf)
% 2 - Ntot_0 = initial number of particles at t=0. [From defineIC: N_tot_0_all.(sp)=InitSize_all(sp)/x2_bar]
N_tot_t_all = N_tot_ICnorm.* P.N_tot_0_all.(sp);
Size_tot = m1_x2_RESC0to1 .* N_tot_t_all;

% - From size to biomass [From definePBEParams: InitSize_all = InitBiom_all./sludgeDensity]
Biomass_tot = Size_tot*P.sludgeDensity; % cellDensity is in g/m^3

if P.tspanCut 
 [m1_x2_RESC0to1,m2_x2_RESC0to1,Marg_N_x1_RESC0to1,Marg_N_x2_RESC0to1,N_tot_t_all,Biomass_tot] = fillIn(m1_x2_RESC0to1,m2_x2_RESC0to1,Marg_N_x1_RESC0to1,Marg_N_x2_RESC0to1,N_tot_t_all,Biomass_tot);
end
return

function aSt = from_N_to_pdf_pvf(aSt,sp)
 % Here Marg_N_x2_RESC0to1 are normalized number of particles, so that the sum(Marg_N_x2_RESC0to1) = 1=zeroth moment
 % So I can get the pdf by dividing by the respective widths
 aSt.(sp).Marg_n_pdf_x2_RESC0to1 = aSt.(sp).Marg_N_x2_RESC0to1./aSt.(sp).x2_widths_TRUE_micronCube(:)';
 aSt.(sp).Marg_n_pvf_x2_RESC0to1 = aSt.(sp).Marg_n_pdf_x2_RESC0to1.*aSt.(sp).x2_all_TRUE_micronCube(:)';
return

function varargout = fillIn(varargin)
% fill in the zeros for the steady state part
 lastIdx = max(find(any(varargin{1},2)));
 for jj = 1:nargin
    aVar = varargin{jj};
    nRows = size(aVar,1)-lastIdx;
    aVar(lastIdx+1:end,:) = repmat(aVar(lastIdx,:),nRows,1);
    varargout{jj} = aVar;
 end
return

function m = calcMoment(N_TI,nChar)
 M_NisTilde = true;
 % - integrate along the rows (x2_coord)
 Marg_x1_TI = sum(N_TI,1);
 % - Integral across the columns (x1 coord)
 m = computeMarg(Marg_x1_TI,M_NisTilde,[],[],nChar,'x2');
return





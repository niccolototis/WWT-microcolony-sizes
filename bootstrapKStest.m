clear all 
close all

% From [https://asaip.psu.edu/articles/beware-the-kolmogorov-smirnov-test/]
% The Kolmogorov-Smirnov (KS) test is used in over 500 refereed papers each
% year in the astronomical literature. It is a nonparametric hypothesis
% test that measures the probability that a chosen univariate dataset is
% drawn from the same parent population as a second dataset (the two-sample
% KS test)

global T_sim plotPVF folderData

P = [];
plotPVF = true;

% time simulation parameters 
T_sim = 50;%days. It should be a multiple of 0.5 = operation cycle length 

% ancillary scripts, order is important
P = GeneralSettings(P); 
P = defineIndexes(P);
P = setOperationParams(P);
% - load ASM1seq data
load(strcat([folderData,'P.mat']),'P');
% update P if I am changing options for plotting like 
P = PBEsettings(P);
plotsVisible = false;
[aSt_options,P] = PBEParams(P,plotsVisible); 
P.runASM1seq = false;
varNames = {'t_res','X_c_res','S_c_res','V_res'};
isBigFile = false(size(varNames));
if strcmp(P.parSetup,'chosenSetASM1seq')
 P = changeParam(P,1,P.parToChange);
 [t_all_ASM1seq,X_c_res,~,V_res] = loadData(varNames,isBigFile,[]);
else 
 [t_all_ASM1seq,X_c_res,~,V_res] = loadData(varNames,isBigFile,'BASE');
end
Xgrams_ASM1seq = X_c_res.*V_res;


% - load PBE data
% for i = 1:length(P.lambdasOpt)
i = 1;

P.thisObjFunType = P.objFunTypes{i+1};
lambdaOpt = P.lambdasOpt(i);
load(strcat([folderData,'aSt_',P.thisObjFunType,'_lambda',num2str(lambdaOpt),'_T_sim',num2str(T_sim),'_nCells',num2str(P.nCells),'.mat']),'aSt');
plotPVF = 1;
aSt = prepareToPlot(aSt,aSt_options);

plot(aSt.AUT.pl.x2_b_TRUE_micronCube,ones(size(aSt.AUT.pl.x2_b_TRUE_micronCube)),'ro')
hold on
plot(aSt.AUT.pl.x2_all_TRUE_micronCube,ones(size(aSt.AUT.pl.x2_all_TRUE_micronCube)),'bx')

x = aSt.AUT.pl.x2_all_TRUE_micronCube;
if plotPVF
 simul = aSt.AUT.pl.Marg_n_pvf_x2_RESC0to1; 
 simul = simul(end,:);
 if strcmp(P.initialPdf,'x2MarginalByKDEstimate') 
    ref_data = aSt.AUT.pl.n_pvf_data(1,:);
    ref_qH = aSt.AUT.pl.n_pvf_qH(1,:);
    ref_qL = aSt.AUT.pl.n_pvf_qL(1,:);
 end
end

%% Generating samples from discrete distribution
addpath('./sampleFromDist')

close all

% here I treat the pdf as height of histograms
plot(1:length(ref_data),ref_data)

ref_data_N = ref_data/sum(ref_data);
simul_N = simul/sum(simul);

% Check that the shape is conserved
% figure
% plot(1:length(ref_data),ref_data_N,'-r')

sizes_sampled_from_ref_data = discretesample(ref_data_N, 1000000);
sizes_sampled_from_simul = discretesample(simul_N, 1000000);

figure
histogram(sizes_sampled_from_ref_data)
hold on
histogram(sizes_sampled_from_simul)

%%

[h,p,ks2stat] = kstest2(sizes_sampled_from_ref_data,sizes_sampled_from_simul,'Alpha',0.05);

% Hypothesis test result, returned as a logical value.
% If h = 1, this indicates the rejection of the null hypothesis at the Alpha significance level.
% If h = 0, this indicates a failure to reject the null hypothesis at the Alpha significance level.
fprintf('Significantly different (null hypot rejected) %s \n',string(h))

% Asymptotic p-value of the test, returned as a scalar value in the range
% (0,1). p is the probability of observing a test statistic as extreme as,
% or more extreme than, the observed value under the null hypothesis. The
% asymptotic p-value becomes very accurate for large sample sizes, and is
% believed to be reasonably accurate for sample sizes n1 and n2, such that
% (n1*n2)/(n1 + n2) ≥ 4.
fprintf('Computed p value p = %g \n',p)

% Test statistic, returned as a nonnegative scalar value.
fprintf('Distance = %g \n',ks2stat)

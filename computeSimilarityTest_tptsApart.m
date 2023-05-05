% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Computation of the three statistics KS AD CVM tests to assess the
% probability that the null hypothesis h0 "The sample of microcolonies I
% got experimentally is drawn from the sequence of distributions generated
% from my model" is true. 
% NOTES: 
% 1 - with *_tptsApart* the statisctic for
% AD CVM tests is saved for each different time point -> one p-value is
% computed for each time point, then an average p-value is considered.
% 2 - To compute the distances and p-values for AD and CVM, a probability
% distribution object for each NDF vector at each time point needs to be
% generated
% 3 - These compPD objects need to be used to generate new samples. Because
% [expSample -> compPD estimation -> randSample] introduces approximations,
% then expSample cannot be really compared to the generated samples ->
% p-values becomes inaccurate. Solution: I am replacing expSample by the
% generated sample with minimum distance from the compPD object at t = 0
% (initial onditions). This is done in section 2 of
% compute_crit_p_tptsApart()
% 4 - Weights depend on the sample. Each element in the sample has its own
% weight that depends on which grid bin does the element belongs to. = >
% every new random sample has a new different weight
% Author: Niccolo Totis
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% With *_tptsAllin*, instead, the distances for all time points
% are summed = > global distance => global_distance_expSample is compared to
% global_distance_MCrandom samplesa

function Dist = computeSimilarityTest_tptsApart(aSt,Lambda_i)
global alpha
close all 

% Small values of p cast doubt on the validity of the null hypothesis.
% Null Hyp: sample comes from the same simulated distribution 
% low p = > sample hardly comes from the same distribution => 2 distrib
% significantly different

alpha = 0.05;

% % save('./computeSimilarityStart.mat')
if isempty(gcp('nocreate'))
 parpool([1 36])
end

%% Choose how is the experimental sample defined. 
% Can be either the original dataset (expSizes), or an artificial sample
% created from the first slice od the NDF. I am following either one of the
% following processes:
% a] 'sampledFrom1stNDF': expSizes -*> expPD -> expCDF -> NDF -> mysampling -*> compPD -> compPDF 
% b] 'expSizes': expSizes -> expPD -> expPDF
howIsExpSampleDefined = 'sampledFrom1stNDF'; % {'sampledFrom1stNDF','expSizes'}
% compPD is evaluated at expSample, which is either a] or b] above


sp = 'AUT';
aSt = defineObjWeights(aSt,sp);
aStSp = aSt.(sp);
x2_ind_toPlot = aStSp.x2_ind_toPlot;
N_res = aStSp.Marg_N_x2_RESC0to1;

nTpoints = size(N_res,1);

x2_all = aStSp.x2_all_TRUE_micronCube(:)';
x2_b = aStSp.x2_b_TRUE_micronCube(:)';
x2_lb = x2_b(1:end-1);
x2_ub = x2_b(2:end);
x2_widths = aStSp.x2_widths_TRUE_micronCube(:)';
myAddedWeights = aStSp.weights;


%% Compute eCDF of the sample
% load sizes
tt = loadSizes();
expSizes = tt.expSizes;


%% AD and CVM are based on CDF -> Check that the following shapes are close:
% checkNDFtoPDF(expSizes,N_res,x2_all,x2_b,x2_lb,x2_ub,x2_widths)


%% From the discrete probability distribution = > continuous distribution that needs to be evaluated 
whatWeights = {'equal','NDF','PDF'};


%% Chose first method
valMCTol = 0.01;
for nW = 1:length(whatWeights)
 aW = whatWeights{nW};
 [mydist,crit,p] = compute_crit_p_tptsApart(expSizes,howIsExpSampleDefined,x2_all,x2_b,x2_lb,x2_ub,x2_widths,N_res,x2_ind_toPlot,nTpoints,myAddedWeights.(aW),aW,valMCTol,Lambda_i);
 Dist.(aW).mydist = mydist;
 Dist.(aW).crit = crit;
 Dist.(aW).p = p;
end 
return

function [dist_expSample,crit_allTps,p_allTps] = compute_crit_p_tptsApart(expSizes,howIsExpSampleDefined,x2_all,x2_b,x2_lb,x2_ub,x2_widths,N_res,ind_toPlot,nTpoints,myWeights_aW,aW,valMCTol,Lambda_i)
% Here differenlty from the previous function I am taking what is in theory
% the right approach for AD and CVM, meaning to have one p value for all
% time points (NDF slices) and not having nTpoins as many p-vals and
% averaging them. The problem with this approach, although correct, is that
% gives me back a 0 p-value for all lambdas. Which is supposedly correct,
% because experimental and computed distribution are too far apart to be
% considered the same but does not help me stratify the lambdas and pick a
% preferred one.

global kernelWidth plotsWhileStatTests alpha nDecimalsToIncludeNDFsampling folderOutput

warning('off','all')

% NB DURING THE ESTIMATION I AM KEEPING 
% Plot of how the NDF evolved through the simulation
if plotsWhileStatTests
 figure
 hold on
 for mm = 1:size(N_res,1)
    if mm ==1
        plot(x2_all,N_res(mm,:),'linewidth',2)
    else
        plot(x2_all,N_res(mm,:))
    end
 end
 title('base NDF')
end

%% 1] for each time point: sample from NDF + Computing compPD object
% prendo N_res -> create its "sample" vector S: if at x1 NDF = 2.33456 and
% nDecimalsToInclude = 2 → I put 233 times the element x1 in S. I do the
% same for all grid points

fold_S_compPD = [folderOutput 'S_compPD/'];
if ~exist(fold_S_compPD,'dir')
 mkdir(fold_S_compPD)
end
similarityFold = [folderOutput 'computeSimilarityTests/'];
if ~exist(similarityFold,'dir')
 mkdir(similarityFold)
end

labLambda = strrep(num2str(Lambda_i),'.','p');
if ~isfile([fold_S_compPD 'S_compPD_lambda_' labLambda '.mat'])
 reestimNDFfig = figure('Visible','on');
 hold on
 compPD = [];
 nPartStart_all = nan(size(N_res(:,:)));

 for oneTpoint = 1:nTpoints
     aTp = ['t' num2str(oneTpoint)];
    oneNDF = N_res(oneTpoint,:);
    nPartStart_all(oneTpoint,:) = round(oneNDF*10^nDecimalsToIncludeNDFsampling);
    nPartThisSlice = sum(nPartStart_all(oneTpoint,:));
    if nPartThisSlice<15 
        % to avoid that I end up with a too small sample (for instance if the distribution is very flat close to zero)
        % NOTE I have already verified that is important that it is
        % important that nDecimalsToIncludeNDFsampling is the same for all
        % time points
        error('Increase nDecimalsToIncludeNDFsampling by one, there for some time points the sampling from the NDF ends up in too little particles to do a proper PD estimation ')
    end

    plotDistributionNDFSamples = false;
    if plotDistributionNDFSamples
        oneSampleFig = figure;
        hold on
    else
        oneSampleFig = [];
    end
    aS = [];

    parfor i = 1:length(oneNDF)
        % take the number of particles to be sampled for each grid interval/cell
        nPartStartThisInt = nPartStart_all(oneTpoint,i);
        if nPartStartThisInt>0
            % random sampling n-2 particles. The last 2 are added to make
            % the meanThisSample on point
            sampledPartThisInt = sampleThisInt(nPartStartThisInt,x2_all,x2_lb,x2_ub,x2_widths,i,oneTpoint,oneSampleFig,plotDistributionNDFSamples);

            % add sampledPartThisInt to the main sample
            aS = [aS sampledPartThisInt];
        end
    end

    S.(aTp) = sort(aS(:));
    [compPD_aTp,N_drawn] = compute_pd_NDF_PDF(S.(aTp),x2_b,[],kernelWidth,{'pd','NDF'});
    % Verify how did the estimaiton go
    figure(reestimNDFfig)
    if oneTpoint ==1
        plot(x2_all,N_drawn,'linewidth',2)
    else
        plot(x2_all,N_drawn)
    end
    compPD.(aTp) = compPD_aTp;
    if oneTpoint ==1 && ~isfile([fold_S_compPD 'S_compPD_t0.mat'])
        compPD_t0 = compPD_aTp;
        S_t0 = S;
        save([fold_S_compPD 'S_compPD_t0.mat'],'S_t0','compPD_t0')
    end
 end
 figure(reestimNDFfig)
 title('sampled reestimated NDF')
 % save
 save([fold_S_compPD 'S_compPD_lambda_' labLambda '.mat'],'S','compPD')
 saveas(reestimNDFfig,[similarityFold 'reestimNDFfig_lambda_' labLambda '.fig'])
else
 load([fold_S_compPD 'S_compPD_lambda_' labLambda '.mat'],'S','compPD')
end

%% Now I resize to just the relevant part (ind_toPlot) and truncate compPD
x2_all = x2_all(ind_toPlot);
x2_lb = x2_lb(ind_toPlot);
x2_ub = x2_ub(ind_toPlot);
x2_b = x2_b([ind_toPlot ind_toPlot(end)+1]);
x2_widths = x2_widths(ind_toPlot);
x2_min = min(x2_b);
x2_max = max(x2_b);

% Truncate the probability distributions!
for oneTpoint = 1:nTpoints
 aTp = ['t' num2str(oneTpoint)];
 compPD.(aTp) = truncate(compPD.(aTp),x2_min,x2_max);
end
 
%% 2] select expSample - my expSample is going to be a fake one
% It is generated by taking the median sample of n = 1000 sample randomly
% generated from the first NDF. This avoids the problem that the estimation
% of the compPD: expSample -> compPD -> randomSample introduces some
% inaccuracies that might cause AD/CVM statistics of randomSamples to not
% be as they should (i.e. evenly 50% 50% smaller or bigger than the
% expSample, given that compPD is generated only from expSample)
switch howIsExpSampleDefined
 case 'expSizes'
    % In this case 
    expSizes = expSizes(expSizes>min(x2_b) & expSizes<max(x2_b));
    expSample = sort(expSizes);
 case 'sampledFrom1stNDF'
    % In this case I am taking the first sample from the NDF as the experimental
    expSample = S.t1;        
end
n = length(expSample);
% Just a Check
[~,hc] = computeWeightsForThisSample(expSample,x2_b,myWeights_aW);
hl = round(N_res(1,ind_toPlot)*10^nDecimalsToIncludeNDFsampling);
if ~isequal(hc,hl)
 disp ([hl' hc'])
 error('Should be the same')
end

tmpFname = [similarityFold 'tmp_expSample_weights_' aW '.mat'];
% if ~isfile(tmpFname)
 n = length(expSample); expSample=[];
 % Creo n = 1000 fake expSamples e poi prendo la moda fra questi
 nF = 2000;% number of fake samples
 [KS_fakeExpSamples,AD_fakeExpSamples,CVM_fakeExpSamples] = deal(zeros(nF,1));
 [randSampleS] = deal(zeros(nF,n));
 [N_drawnS] = deal(zeros(nF,length(x2_all)));
 parfor rep = 1:nF
    if mod(rep,100) ==0
        fprintf('Lambda %g : fake expSample %g \n',Lambda_i, rep)
    end

    % draw n random samples from the compPD objects
    randSample = [];
    m = n;
    while m>0
        randSampleTmp = random(compPD.t1,m,1);
        dd = find(randSampleTmp<x2_min & randSampleTmp>x2_max);
        if ~isempty(dd)
            error('WTF I have truncated the pd')
        end
        randSampleTmp = randSampleTmp(randSampleTmp>x2_min & randSampleTmp<x2_max);
        randSample = [randSample; randSampleTmp];
        m = n-length(randSample);
    end
    randSample = sort(randSample);
    randSampleS(rep,:) = randSample;
    [~,N_drawn] = compute_pd_NDF_PDF(randSample,x2_b,[],0.5,{'pd','NDF'});
    N_drawnS(rep,:) = N_drawn;
    compCDF_at_randSample = cdf(compPD.t1,randSample);
    [myAddedWeights_aW_,~] = computeWeightsForThisSample(randSample,x2_b,myWeights_aW);

    % Store
    [KS_fakeExpSample_rep,~] = computeStat('KS',compCDF_at_randSample,n,myAddedWeights_aW_,randSample,compPD.t1);
    KS_fakeExpSamples(rep) = KS_fakeExpSample_rep;
    AD_fakeExpSamples(rep) = computeStat('AD',compCDF_at_randSample,n,myAddedWeights_aW_);
    CVM_fakeExpSamples(rep) = computeStat('CVM',compCDF_at_randSample,n,myAddedWeights_aW_);
 end      
 % The randomly generated sample which has the minimum distance from the
 % distribution it has been generated from is the one which best
 % represents the distribution itself. Thus, I take this one as my
 % expSample. As a consequence 99.99% percent of other random samples have bigger
 % distances p-val = 1 that Ho (expSample comes from the distribution) is true
 [~,whereKS] = min(KS_fakeExpSamples);
 [~,whereAD] = min(AD_fakeExpSamples);
 [~,whereCVM] = min(CVM_fakeExpSamples);

 expSample.KS = randSampleS(whereKS,:);
 expSample.AD = randSampleS(whereAD,:);
 expSample.CVM = randSampleS(whereCVM,:);

 save(tmpFname)
% else
% load(tmpFname,'expSample','nF','N_drawnS','whereKS','whereAD','whereCVM')
% end
% recap plot
figure
hold on
for rep = 1:nF
 plot(x2_all,N_drawnS(rep,:))
end
plot(x2_all,N_drawnS(whereKS,:),'g','linewidth',2)
plot(x2_all,N_drawnS(whereAD,:),'r','linewidth',2)
plot(x2_all,N_drawnS(whereCVM,:),'b','linewidth',2)
plot(x2_all,N_res(1,ind_toPlot),'--k','linewidth',1)


%% 3] Compute the weights for the expSample
[myAddedWeights_aW_x,~] = computeWeightsForThisSample(expSample.KS,x2_b,myWeights_aW);
myAddedWeights_aW.KS = myAddedWeights_aW_x;
[myAddedWeights_aW_x,~] = computeWeightsForThisSample(expSample.AD,x2_b,myWeights_aW);
myAddedWeights_aW.AD = myAddedWeights_aW_x;
[myAddedWeights_aW_x,~] = computeWeightsForThisSample(expSample.CVM,x2_b,myWeights_aW);
myAddedWeights_aW.CVM = myAddedWeights_aW_x;

%% 4] For all time points use compPD objects to compute compCDF at expSample, then compute the 3 statistics
% For KS test, directly compute the p-value
figure

hold on
plot(expSample.KS,0.05,'g*')
plot(expSample.AD,0.05,'rx')
plot(expSample.CVM,0.05,'bo')

% statisctic for KS test can only be calculated Allin = > no need to store
% all the different time points
% p.KS gives me the probability that the null is REJECTED. Close to zero
% the null is true
dist_expSample.KS = 0; p.KS=0; crit.KS=NaN;

% statisctic for AD CVM tests is saved for each different time point -> one
% p-value is computed for each time point, then an average p-value is
% considered
dist_expSample.KS = 0; 
dist_expSample.AD = zeros(nTpoints,1); 
dist_expSample.CVM = zeros(nTpoints,1);

compCDF_at_expSample_allTpts.KS = nan(length(expSample.KS),nTpoints);
compCDF_at_expSample_allTpts.AD = nan(length(expSample.AD),nTpoints);
compCDF_at_expSample_allTpts.CVM = nan(length(expSample.CVM),nTpoints);

figure
hold on
for oneTpoint = 1:nTpoints
 aTp = ['t' num2str(oneTpoint)];
 
 compCDF_at_expSample_allTpts.KS(:,oneTpoint) = cdf(compPD.(aTp),expSample.KS);   
 compCDF_at_expSample_allTpts.AD(:,oneTpoint) = cdf(compPD.(aTp),expSample.AD);   
 compCDF_at_expSample_allTpts.CVM(:,oneTpoint) = cdf(compPD.(aTp),expSample.CVM);     
 
 % @plots

 plot(expSample.KS,compCDF_at_expSample_allTpts.KS(:,oneTpoint),'g-','linewidth',1)
 plot(expSample.AD,compCDF_at_expSample_allTpts.AD(:,oneTpoint),'r-','linewidth',1)
 plot(expSample.CVM,compCDF_at_expSample_allTpts.CVM(:,oneTpoint),'b-','linewidth',1)

 % KS
 [dist_expSample_KS_aTp,pKS_aTp] = computeStat('KS',compCDF_at_expSample_allTpts.KS(:,oneTpoint),n,myAddedWeights_aW.KS,expSample.KS,compPD.(aTp));
 % KS is the only statistics that gives me the p-value without the need
 % to perform a MC approach = > I can compute the p-value for all time
 % points immediately
 dist_expSample.KS = dist_expSample.KS+dist_expSample_KS_aTp;
 p.KS = p.KS+pKS_aTp;
 
 % AD
 dist_expSample.AD(oneTpoint) = computeStat('AD',compCDF_at_expSample_allTpts.AD(:,oneTpoint),n,myAddedWeights_aW.AD);
 
 % CVM
 dist_expSample.CVM(oneTpoint) = computeStat('CVM',compCDF_at_expSample_allTpts.CVM(:,oneTpoint),n,myAddedWeights_aW.CVM);
end
hold off

% ONLY for KS test the distance and the p values are averaged (already
% summed thoruhout the loop)
dist_expSample.KS = dist_expSample.KS/nTpoints;
p.KS = p.KS/nTpoints;
crit.KS = NaN;

%% 5] Monte Carlo computation of p values and critical values: ONLY for AD and CVM 
% Parameters for MC determination of p-values
mctol = valMCTol; % Standard error
vartol = mctol^2; % standard desired variance

for oneTpoint = 1:nTpoints
 aTp = ['t' num2str(oneTpoint)];
 
 % Empty distance for all time points, p-val and crit-val
 AD_MC = []; CVM_MC=[];
 crit.AD.(aTp) = 0; p.AD.(aTp)=0;
 crit.CVM.(aTp) = 0; p.CVM.(aTp)=0;
 
 mcRepsTot = 0; % Counting the MC reps performed so far for both AD and CVM
 mcRepsMin = 1000; % Min MC reps to be performed for AD and CVM
 while true
    mcRepsOld = mcRepsTot; % MC performed so far for both AD and CVM
    mcReps = ceil(mcRepsMin - mcRepsOld); % remaining MC reps to be performed
    mcRepsTot = mcRepsOld + mcReps;

    % Compute aStatMC
    parfor rep = 1:mcReps
        if mod(rep,100) ==0
            fprintf('Lambda %g Weight %s : tPoint %g MCrep = %g \n',Lambda_i,aW,oneTpoint, rep)
        end

        % draw n random samples from the compPD objects
        randSample = [];
        m = n;
        while m>0
            randSampleTmp = random(compPD.(aTp),m,1);
            dd = find(randSampleTmp<x2_min & randSampleTmp>x2_max);
            if ~isempty(dd)
                error('WTF I have truncated the pd')
            end
            randSampleTmp = randSampleTmp(randSampleTmp>x2_min & randSampleTmp<x2_max);
            randSample = [randSample; randSampleTmp];
            m = n-length(randSample);
        end
        randSample = sort(randSample);
        compCDF_at_randSample = cdf(compPD.(aTp),randSample);
        [myAddedWeights_aW,~] = computeWeightsForThisSample(randSample,x2_b,myWeights_aW);
        
        % Store
        AD_MC(oneTpoint,rep) = computeStat('AD',compCDF_at_randSample,n,myAddedWeights_aW);
        CVM_MC(oneTpoint,rep) = computeStat('CVM',compCDF_at_randSample,n,myAddedWeights_aW);
    end
    
    critMC.AD.(aTp) = prctile(AD_MC(oneTpoint,(end-mcReps+1):end), 100*(1-alpha));
    pMC.AD.(aTp) = sum(AD_MC(oneTpoint,(end-mcReps+1):end) > dist_expSample.AD(oneTpoint))./mcReps;
    crit.AD.(aTp) = (mcRepsOld*crit.AD.(aTp) + mcReps*critMC.AD.(aTp)) / mcRepsTot;
    p.AD.(aTp) = (mcRepsOld*p.AD.(aTp) + mcReps*pMC.AD.(aTp)) / mcRepsTot;
    sepsq.AD.(aTp) = max(p.AD.(aTp)*(1-p.AD.(aTp))/mcRepsTot, 1/mcRepsTot^2);

    critMC.CVM.(aTp) = prctile(CVM_MC(oneTpoint,(end-mcReps+1):end), 100*(1-alpha));
    pMC.CVM.(aTp) = sum(CVM_MC(oneTpoint,(end-mcReps+1):end) > dist_expSample.CVM(oneTpoint))./mcReps;
    crit.CVM.(aTp) = (mcRepsOld*crit.CVM.(aTp) + mcReps*critMC.CVM.(aTp)) / mcRepsTot;
    p.CVM.(aTp) = (mcRepsOld*p.CVM.(aTp) + mcReps*pMC.CVM.(aTp)) / mcRepsTot;
    sepsq.CVM.(aTp) = max(p.CVM.(aTp)*(1-p.CVM.(aTp))/mcRepsTot, 1/mcRepsTot^2);

    sepsq_both = max(sepsq.AD.(aTp),sepsq.CVM.(aTp));
    if sepsq_both < vartol
        break
    end
    % Based on the current estimate, find the number of trials needed to
    % make the MC std err less than the specified tolerance.
    mcRepsMin = 1.2 * (mcRepsTot*sepsq_both)/vartol;
 end
end

%% media
p_allTps.KS = p.KS;
p_allTps.AD = 0; crit_allTps.AD=0;
p_allTps.CVM = 0; crit_allTps.CVM=0;
for oneTpoint = 1:nTpoints
 aTp = ['t' num2str(oneTpoint)];
 
 % AD
 p_allTps.AD = p_allTps.AD+p.AD.(aTp);
 crit_allTps.AD = crit_allTps.AD+crit.AD.(aTp);
 
 % CVM
 p_allTps.CVM = p_allTps.CVM+p.CVM.(aTp);
 crit_allTps.CVM = crit_allTps.CVM+crit.CVM.(aTp);
end
% Average
p_allTps.AD = p_allTps.AD/nTpoints;
crit_allTps.AD = crit_allTps.AD/nTpoints;

p_allTps.CVM = p_allTps.CVM/nTpoints;
crit_allTps.CVM = crit_allTps.CVM/nTpoints;
return

function expSizes = loadSizes()
filenamesave = strcat(['./outputs/data/MiCol_sizes.mat']);
if isfile(filenamesave)
 expSizes = load(filenamesave);
else
 str = 'untreated_AOB_NOB_NT';
 data = readtable('../exp_data_miCol/200207_200228_FISH_Izico_untreated_treated_Results.xlsx','Sheet',str);
 areas = table2array(data(:,{'Area'}));
 % - remove N/A from data
 areas = areas(~isnan(areas));
 % - sort the data. areas.(str) is a matlab stucture, like a wardrobe with different
 % - shelves, one for each (str)
 areas = sort(areas);
 diams = 2*sqrt(areas./pi); % from areas to diameters
 expSizes = areaToVol(areas);
 filenamesave = strcat(['./outputs/data/MiCol_sizes.mat']);
 save(filenamesave,'expSizes');
end
return

function varargout = computeStat(whatStat,compCDF_at_sample,n,myAddedWeights_aW,varargin)
 global alpha
 myAddedWeights_aW = myAddedWeights_aW(:);
 switch whatStat
    case 'KS'
        expSample = varargin{1};
        compPD_aTp = varargin{2};
        tail = 'unequal'; % both directions
        [h_aTp1,pKS_aTp1,statKS_aTp1,cv_aTp1] = kstest(expSample,'CDF',compPD_aTp,'tail',tail,'alpha',alpha);
        [h_aTp2,pKS_aTp2,statKS_aTp2,cv_aTp2] = kstest(expSample,'CDF',[expSample(:) compCDF_at_sample(:)],'tail',tail,'alpha',alpha);
        if ~isequal(pKS_aTp1,pKS_aTp2) || ~isequal(statKS_aTp1,statKS_aTp2)
            error('both inputs to kstest should be fine')
        end
        aStat = statKS_aTp1;
        varargout{2} = pKS_aTp1;
    case 'AD'
        compCDF_at_expSample = compCDF_at_sample;
        compCDF_at_expSample(compCDF_at_expSample ==1)=0.99999;
        w = 2*(1:n) - 1 ;
        w = w.*myAddedWeights_aW(:)';
        numer = (log(compCDF_at_expSample)+ log(1-compCDF_at_expSample(end:-1:1)));
        aStat = -w*numer/n - n;
        if isinf(aStat)
            pausa_qui = 3;
        end
    case 'CVM'
        compCDF_at_expSample = compCDF_at_sample;
        w = (2*(1:n)' - 1) ./ (2*n);
        w = w.*myAddedWeights_aW(:);
        aStat = 1/(12*n) + sum((w - compCDF_at_expSample).^2);
        if isinf(aStat)
            pausa_qui = 3;
        end
 end
 varargout{1} = aStat;
return

function aSt = defineObjWeights(aSt,sp)
% - Defining the weights that weigh the difference data-simResult in the
% obj function: (yOneParam-yData).*objWeights. Defined inversely
% proportional to the measure of uncertainty (widht of the confidence interval)
 aSt.(sp).weights.equal = ones(size(aSt.(sp).N_qH));
 aSt.(sp).weights.NDF = 1./(aSt.(sp).N_qH-aSt.(sp).N_qL);
 aSt.(sp).weights.PDF = 1./(aSt.(sp).n_pdf_qH-aSt.(sp).n_pdf_qL);
return

function s = areaToVol(a)
 s = ((4/3).*pi.*(sqrt(a/pi)).^3);
return

function sampledPartThisInt = sampleThisInt(nPartStartThisInt,x2_all,x2_lb,x2_ub,x2_widths,i,oneTpoint,oneSampleFig,plotDistributionNDFSamples)
 gridPoint = x2_all(i); sigma=x2_widths(i)/4;
 x2_lb_i = x2_lb(i); x2_ub_i=x2_ub(i);
 
 sampledPartThisInt = NaN;
 
 % SIMPLEST
% sampledPartThisInt = ones(1, nPartStartThisInt)*gridPopint;
 
 % RANDOM SAMPLE inside the interval
 % take the number of particles to be sampled for each grid interval/cell
 if nPartStartThisInt>0
    nPartsWellPut = min(1,nPartStartThisInt);
    % random sampling n-2 particles. The last 2 are added to make
    % the meanThisSample on point
    nThisSample = 0;  unifSampleThisInt=[];
    nPartsLeft = (nPartStartThisInt-nPartsWellPut)-nThisSample;

    while nPartsLeft>0 
        unifSampleThisInt = [unifSampleThisInt normrnd(gridPoint,sigma,1,nPartsLeft)];
        sampledPartThisInt = unifSampleThisInt(unifSampleThisInt<x2_ub(i) & unifSampleThisInt>x2_lb(i));
        sampledPartThisInt = sort(sampledPartThisInt);
        nThisSample = length(sampledPartThisInt);
        nPartsLeft = (nPartStartThisInt-nPartsWellPut)-nThisSample;

        % In case I already put too many particles I remove them,
        % so that I then well put the last particle. Here I am
        % removing the particles furthest away from the center
        if nPartsLeft < 0
            nToRem = abs(nPartsLeft);
            diffs = abs(sampledPartThisInt-gridPoint);
            [~,posx] = maxk(diffs,nToRem);
            sampledPartThisInt(posx) = [];
        end
        nThisSample = length(sampledPartThisInt);
        nPartsLeft = (nPartStartThisInt-nPartsWellPut)-nThisSample;
        meanThisSample = mean(sampledPartThisInt); 
    end
    
    if  plotDistributionNDFSamples
        figure(oneSampleFig)
        plot(sampledPartThisInt,i,'rx')
        xline(x2_lb_i)
        xline(x2_ub_i)
        xline(gridPoint)
    end
    
    sampledPartThisInt(isnan(sampledPartThisInt)) = [];
    
    culo = [0 0 0];
    for nn = 1:nPartsWellPut
        % I add n new particles as long as the grid point is not
        % representative of the whole interval
        % I set the size of the new particle so that the mean in the
        % interval matches the grid point
        newPart = gridPoint*(nThisSample+1)-sum(sampledPartThisInt);
        % If newPart falls out of the boundaries I make it so that it lays
        % inbetween (not exactly on the boundary, just in case)
        if newPart<x2_lb(i)
            newPart = x2_lb(i)+x2_widths(i)/1000;
            if nn ==nPartsWellPut 
               warning(['The sampled particles in tPoint ' num2str(oneTpoint) ' interval ' num2str(i) ' have a mean different from the grid point'])
            end
        end
        if newPart>x2_ub(i)
            newPart = x2_ub(i)-x2_widths(i)/1000;
            if nn ==nPartsWellPut 
                warning(['The sampled particles in tPoint ' num2str(oneTpoint) ' interval ' num2str(i) ' have a mean different from the grid point'])
            end
        end
        culo(nn) = 1;
        plot(newPart,i,'o','MarkerFaceColor',culo)
        sampledPartThisInt = [sampledPartThisInt newPart];
        nThisSample = nThisSample+1;
    end       
 end
return

function [myAddedWeights_aW,hc] = computeWeightsForThisSample(sample,x2_b,myWeights_aW)
 % As experimental CDF (eCDF) I do not need to compute it by myself I only
 % need to evaluate with cdf() the PDobj of the first slice
 hc = histcounts(sample,x2_b);
 myAddedWeights_aW = [];
 for aBin = 1:length(hc)
    myAddedWeights_aW = [myAddedWeights_aW repmat(myWeights_aW(aBin),1,hc(aBin))];
 end
return

% function [dist_expSample,crit,p] = compute_aStat_critPvals(statsToCompute,expSizes,howIsExpSampleDefined,x2_all,x2_b,x2_lb,x2_ub,x2_widths,N_res,ind_toPlot,nTpoints,myWeights_aW,valMCTol,Lambda_i)
% global kernelWidth plotsWhileStatTests alpha nDecimalsToIncludeNDFsampling
% x2_all = x2_all(ind_toPlot);
% x2_lb = x2_lb(ind_toPlot);
% x2_ub = x2_ub(ind_toPlot);
% x2_b = x2_b([ind_toPlot ind_toPlot(end)+1]);
% x2_widths = x2_widths(ind_toPlot);
% x2_min = min(x2_b);
% x2_max = max(x2_b);
% 
% % Plot of how the NDF evolved through the simulation
% if plotsWhileStatTests
% figure
% hold on
% for mm = 1:size(N_res,1)
%     if mm ==1
%         plot(x2_all,N_res(mm,ind_toPlot),'linewidth',2)
%     else
%         plot(x2_all,N_res(mm,ind_toPlot))
%     end
% end
% end
% 
% %% 1] for each time point: sample from NDF + Computing compPD object
% % prendo N_res -> create its "sample" vector S: if at x1 NDF = 2.33456 and
% % nDecimalsToInclude = 2 → I put 233 times the element x1 in S. I do the
% % same for all grid points
% for oneTpoint = 1:nTpoints
% aTp = ['t' num2str(oneTpoint)];
% oneNDF = N_res(oneTpoint,ind_toPlot);
% nPartStart_all.(aTp) = round(oneNDF*10^nDecimalsToIncludeNDFsampling);
% nPartThisSlice = sum(nPartStart_all.(aTp));
% if nPartThisSlice<15 
%     % to avoid that I end up with a too small sample (for instance if the distribution is very flat close to zero)
%     % NOTE I have already verified that is important that it is
%     % important that nDecimalsToIncludeNDFsampling is the same for all
%     % time points
%     error('Increase nDecimalsToIncludeNDFsampling by one, there for some time points the sampling from the NDF ends up in too little particles to do a proper PD estimation ')
% end
% 
% aS = [];
% for i = 1:length(oneNDF)
%     % take the number of particles to be sampled for each grid interval/cell
%     nPartStartThisInt = nPartStart_all.(aTp)(i);
%     if nPartStartThisInt>0
%         % random sampling n-2 particles. The last 2 are added to make
%         % the meanThisSample on point
%         unifSampleThisInt = x2_lb(i) + x2_widths(i)*rand(1,nPartStartThisInt-2);
%         sampleThisInt = unifSampleThisInt;
%         nThisSample = length(sampleThisInt);
%         meanThisSample = mean(sampleThisInt);
%         gridPoint = x2_all(i);
%         % I add n new particles as long as the grid point is not
%         % representative of the whole interval
%         nP = 1;
%         while nThisSample<nPartStartThisInt 
%             % fprintf('n = %g particle added to int %g \n',nP,i)
%             % adjust the size of the new particle so that the mean in the
%             % interval matches
%             newPart = gridPoint*(nThisSample+1)-sum(sampledPartThisInt);
%             % If newPart falls out of the boundaries I make it so that it lays
%             % inbetween (not exactly on the boundary, just in case)
%             if newPart<x2_lb(i)
%                 newPart = x2_lb(i)+x2_widths(i)/100;
%             end
%             if newPart>x2_ub(i)
%                 newPart = x2_ub(i)-x2_widths(i)/100;
%             end
%             sampledPartThisInt = [sampledPartThisInt newPart];
%             nThisSample = length(sampledPartThisInt);
%             meanThisSample = mean(sampledPartThisInt);
%             nP = nP+1;
%         end
%         % add sampledPartThisInt to the main sample
%         aS = [aS sampledPartThisInt];
%     end
% end
% S.(aTp) = sort(aS(:));
% compPD_aTp = compute_pd_NDF_PDF(S.(aTp),[],[],kernelWidth,{'pd'});
% compPD.(aTp) = compPD_aTp;
% end
% 
% 
% %% 2] select expSample 
% switch howIsExpSampleDefined
% case 'expSizes'
%     % In this case 
%     expSizes = expSizes(expSizes>min(x2_b) & expSizes<max(x2_b));
%     expSample = sort(expSizes);
% case 'sampledFrom1stNDF'
%     % In this case I am taking as the experimental
%     aTp = ['t' num2str(1)];
%     expSample = S.(aTp);
%     % As experimental CDF (eCDF) I do not need to compute it by myself I only
%     % need to evaluate with cdf() the PDobj of the first slice
%             hc = histcounts(expSample,x2_b);
%             myAddedWeights_aW = [];
%             for aBin = 1:length(x2_all)
%                 myAddedWeights_aW = [myAddedWeights_aW repmat(myWeights_aW(aBin),1,hc(aBin))];
%             end
% end
% n = length(expSample);
% 
% %% 3] Compute compCDF at expSample + the statistics for all time points 
% figure(1)
% 
% hold on
% plot(expSample,0.05,'rx')
% dist_expSample.KS = 0; p.KS=0;
% dist_expSample.AD = 0; 
% dist_expSample.CVM = 0; 
% 
% 
% for oneTpoint = 1:nTpoints
% % IMPO here I compare 
% aTp = ['t' num2str(oneTpoint)];
% compCDF_at_expSample_allTpts(:,oneTpoint) = cdf(compPD.(aTp),expSample);     
% 
% % @plots
% plot(expSample,compCDF_at_expSample_allTpts(:,oneTpoint),'b-','linewidth',1)
% 
% for nStat = 1:length(statsToCompute)
%     whatStat = statsToCompute{nStat};
%     switch whatStat
%         case 'KS'
%             [mydist_aStat_aTp,pKS_aTp] = computeStat(whatStat,compCDF_at_expSample_allTpts(:,oneTpoint),n,myAddedWeights_aW,expSample,compPD.(aTp),compCDF_at_expSample_allTpts(:,oneTpoint));
%             p.KS = p.KS+pKS_aTp;
%         otherwise
%             % Compute the statistic for the experimental sample
%             mydist_aStat_aTp = computeStat(whatStat,compCDF_at_expSample_allTpts(:,oneTpoint),n,myAddedWeights_aW);
%     end
%     mydist.(whatStat).(aTp) = mydist_aStat_aTp;
%     dist_expSample.(whatStat) = dist_expSample.(whatStat)+mydist_aStat_aTp;
% end
% end
% % ONLY for KS test the distance and the p values are averaged (already
% % summed thoruhout the loop)
% dist_expSample.KS = dist_expSample.KS/nTpoints;
% p.KS = p.KS/nTpoints;
% crit.KS = NaN;
% 
% 
% %% 3] Monte Carlo computation of p values and critical values: ONLY for AD and CVM 
% % Compute the same statistic for random samples drawn form the same compPD object
% statsToCompute = setdiff(statsToCompute,'KS','stable');
% 
% p_allTps.AD = 0; p_allTps.CVM=0;
% crit_allTps.AD = 0; crit_allTps.CVM=0;
% for oneTpoint = 1:nTpoints
% aTp = ['t' num2str(oneTpoint)];
% 
% mctol = valMCTol; % Standard error
% vartol = mctol^2; % standard desired variance
% crit_aTp.AD.(aTp) = 0; crit_aTp.CVM.(aTp) = 0;
% p_aTp.AD.(aTp) = 0; p_aTp.CVM.(aTp) = 0;
% mcRepsTot = 0;
% mcRepsMin = 1000;    
% 
% while true
%     mcRepsOld = mcRepsTot;
%     mcReps = ceil(mcRepsMin - mcRepsOld);
%     AD_MC_aTp = zeros(mcReps,1);
%     CVM_MC_aTp = zeros(mcReps,1);
% 
%     % Compute aStatMC
%     parfor rep = 1:mcReps
% %     for rep = 1:mcReps
%         if mod(rep,20) ==0
%             fprintf('Lambda %g : tPt = %g  MCrep=%g \n',Lambda_i, oneTpoint, rep)
%         end
%         xRand = [];
%         m = n;
%         while m>0
%             xRandTent = random(compPD.(aTp),m,1);
%             xRandTent = xRandTent(xRandTent>x2_min & xRandTent<x2_max);
%             xRand = [xRand; xRandTent];
%             m = n-length(xRand);
%         end
%         xRand = sort(xRand);
% 
%         compCDF_at_xRand = cdf(compPD.(aTp),xRand);
% 
%         AD_MC_aTp_aRep = computeStat('AD',compCDF_at_xRand,n,myAddedWeights_aW);
%         CVM_MC_aTp_aRep = computeStat('CVM',compCDF_at_xRand,n,myAddedWeights_aW);
% 
%         AD_MC_aTp(rep) = AD_MC_aTp_aRep;
%         CVM_MC_aTp(rep) = CVM_MC_aTp_aRep;
%     end
%     distMC_aTp.AD = AD_MC_aTp;
%     distMC_aTp.CVM = CVM_MC_aTp;
% 
%     for nStat = 1:length(statsToCompute)
%         whatStat = statsToCompute{nStat};
%         critMC.(whatStat) = prctile(distMC_aTp.(whatStat), 100*(1-alpha));
%         pMC.(whatStat) = sum(distMC_aTp.(whatStat) > mydist.(whatStat).(aTp))./mcReps;
% 
%         mcRepsTot = mcRepsOld + mcReps;
%         crit_aTp.(whatStat).(aTp) = (mcRepsOld*crit_aTp.(whatStat).(aTp) + mcReps*critMC.(whatStat)) / mcRepsTot;
%         p_aTp.(whatStat).(aTp) = (mcRepsOld*p_aTp.(whatStat).(aTp) + mcReps*pMC.(whatStat)) / mcRepsTot;
% 
%         % Compute a std err for p, with lower bound (1/N)*(1-1/N)/N when p ==0.
%         sepsq.(whatStat) = max(p_aTp.(whatStat).(aTp)*(1-p_aTp.(whatStat).(aTp))/mcRepsTot, 1/mcRepsTot^2);
%     end
% 
%     % look for the maximum of sepsq
%     sepsq = -Inf; 
%     remove = [];
%     for nStat = 1:length(statsToCompute)
%         whatStat = statsToCompute{nStat};
%         sepsq = max(sepsq,sepsq.(whatStat));
%         % stop the computation for this specific stat
%         if sepsq.(whatStat) < vartol
%             remove = [remove nStat];
%         end
%     end
%     statsToCompute(remove) = [];
%     if isempty(statsToCompute)
%         break
%     end
%     % Based on the current estimate, find the number of trials needed to
%     % make the MC std err less than the specified tolerance.
%     mcRepsMin = 1.2 * (mcRepsTot*sepsq)/vartol;
% end
% p_allTps.AD = p_allTps.AD+p_aTp.(whatStat).(aTp);
% p_allTps.CVM = p_allTps.CVM+p_aTp.CVM.(aTp);
% crit_allTps.AD = crit_allTps.AD+crit_aTp.AD.(aTp);
% crit_allTps.CVM = crit_allTps.CVM+crit_aTp.CVM.(aTp);
% end
% p.AD = p_allTps.AD/nTpoints;
% p.CVM = p_allTps.CVM/nTpoints;
% crit.AD = crit_allTps.AD/nTpoints;
% crit.CVM = crit_allTps.CVM/nTpoints;
% return


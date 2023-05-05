function [aSt,P] = kdestimation(aSt,P)
global speciesPBE spNames x2_ind_toPlot hotToComputeNDF plotPDFfckernewidth wtbwidthftkerneloexpPD kernelWidth excludeFirstPoint checkECDF biom_lower_than_SS nCells oneGridPointPerSCell sC nTrapezPerCell logBasis x2MinIsSCell cutAt_x2Max logBasis
global wtfresamplearefitdist cutInitialNDFafterEstim cutInitialNDFbeforeEstim IC_at 
% In this script I am estimating N, which is the vector of zeroth moments,
% ie the number of particles allocated in each bin of the coordinate grid.
% This is the NDF integrated over the width of the interval
% Then I am also estimating the PDF, used later on as a weighting factor 

for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE}); 
aSt.(sp).diamSCell = 1;% micron
areaSC = pi*(aSt.(sp).diamSCell/2)^2;% micron
x2SCell = areaToVol(areaSC); % from areas to expSizes=volumes

filenamesave = strcat(['./outputs/data/MiCol_sizes.mat']);
if isfile(filenamesave)
    load(filenamesave);
else
    data = readtable('./exp_data_miCol/AOB_NOB_area_volume.xls');
    expSizes = areasToExpSizes(data);
    save(filenamesave,'expSizes')
end

% I put half of the greed points of nCells between the minimum and maximum measured expSizes
if nCells<4
 ME = MException('Increase nCells');
 throw(ME);
end
switch sC
 case 'lin'
    nCells_half = round(nCells/2);
    x2_all = linspace(x2SCell,max(expSizes),nCells_half)';
    x2_widths = diff(x2_all);
    oneWidth = x2_widths(1);
    x2_all = (x2SCell:oneWidth:max(expSizes)+(oneWidth/2))';
    
    % Add the second half of grid points after sMax
    x2_beforeSC = sort(x2SCell:-oneWidth:min(expSizes));
    x2_beforeSC(x2_beforeSC ==x2SCell)=[];
    nCells_2half = nCells-nCells_half;
    x2_afterMax = (oneWidth:oneWidth:(oneWidth*nCells_2half))';
    x2_all = [x2_beforeSC; x2_all; max(expSizes)+x2_afterMax];
    x2_lb = x2_all-oneWidth/2;
    x2_ub = x2_all+oneWidth/2;
 case 'log'
    expMin = log(min(expSizes))/log(logBasis);
    expSc = log(x2SCell)/log(logBasis);
    expMax = log(3*max(expSizes))/log(logBasis);
    % First i consider exp_b actually at the midpoint,hen I shift it
    % down + I add the last point. In this way I have that expSc and
    % expMax are at the center of expCells
    exp_b = linspace(expSc,expMax,nCells)';
    expCellWidth = mode(diff(exp_b));
    % adding the part between expMin and expSc, spaced as the part
    % between expSc and expMax.I go backwards because I might not end
    % up exatly at the endpoint
    exp_b_firstPart = sort(expSc:-expCellWidth:expMin);
    exp_b = [exp_b_firstPart(1:end-1)'; exp_b];% avoid taking expSc twice
    expCellWidth = mode(diff(exp_b));
    % then I shift it down + I add the last point
    exp_b = exp_b-expCellWidth/2;
    exp_b = [exp_b ; exp_b(end)+expCellWidth];
    exp_lb = exp_b(1:end-1);
    exp_ub = exp_b(2:end);
    exp_all = (exp_lb+exp_ub)/2;
    reallyEqual(exp_all,exp_lb+expCellWidth/2,' = ','Should be the same why they are not')
    x2_all = logBasis.^exp_all;
    x2_b = logBasis.^exp_b;
    x2_lb = logBasis.^exp_lb;
    x2_ub = logBasis.^exp_ub;
end

if x2MinIsSCell
 x2_lb(1) = x2_all(1); % first lower bound coincides with the single cell 
else
 x2_lb(1) = 0; 
end

x2_b = [x2_lb; x2_ub(end)];
x2_widths = diff(x2_b);


if cutInitialNDFbeforeEstim
 [x2_all,x2_b,x2_lb,x2_ub,x2_widths,x2Min,x2Max] = cutDistr(x2_all,x2SCell,x2_b,x2_lb,x2_ub,x2_widths);
end

%% If biom_lower_than_SS
% If this option is active I still do the kdestimation because are the data
% that I use as reference in the movie plot
if strcmp(IC_at,'startupScenario')
 aSt.(sp).mu_x2 = 0; % In this case they are all equal
 aSt.(sp).sigma_x2 = x2SCell/20;
 % - create mu 2x1 and sigma 2x2 matrices later used to build the
 % probability density function in defineIC
 aSt.(sp).mu = [aSt.(sp).mu_x1 aSt.(sp).mu_x2]; % in cm
 aSt.(sp).sigma = [aSt.(sp).sigma_x1 aSt.(sp).sigma_x2].*eye(2); % in cm

 % pdf
 figure
 help = 1000;
 pdf_x2_all = lognpdf(x2_all/x2Max*help, aSt.(sp).mu_x2/x2Max*help, aSt.(sp).sigma_x2*help);
 plot(x2_all,pdf_x2_all);
 title('Initial PDF for low biomass before normalization (kdestimation.m)')
 
 % cdf
 figure
 help = 1000;
 cdf_x2_b = logncdf(x2_b/x2Max*help, aSt.(sp).mu_x2/x2Max*help, aSt.(sp).sigma_x2*help);
 yvals_all = cdf_x2_b(2:end)-cdf_x2_b(1:(end-1));
 plot(x2_all,yvals_all);
 title('Initial NDF for low biomass before normalization (kdestimation.m)')
else
 aSt.(sp).mu_x2 = [];
 aSt.(sp).sigma_x2 = [];
end

hold off
switch hotToComputeNDF
 case 'fitdist'
 %% Estimating PD obj -> PDF -> N using FITDIST
    if wtbwidthftkerneloexpPD 
        % whats the best width for the kernel of the expPD?
        widths = 0.5:0.5:4;
        red = [1, 0, 0];
        green = [0, 1, 0];
        colors_p = [linspace(red(1),green(1),length(widths))', linspace(red(2),green(2),length(widths))', linspace(red(3),green(3),length(widths))'];
        hold on
        for ii = 1:length(widths)
            aWi = widths(ii);
            N_data = compute_pd_NDF_PDF(expSizes,x2_b,[],aWi,{'NDF'});
            plot(x2_all,N_data/max(N_data),'color',colors_p(ii,:));%,'linewidth',1.5)
        end
        legend(cellstr(string(widths)));
    end

    % a] Estimating NDF (N_data) [initial distribution for the simulations]
    [expPD,N_data,PDF_data] = compute_pd_NDF_PDF(expSizes,x2_b,x2_all,kernelWidth,{'pd','NDF','PDF'});

    if plotPDFfckernewidth
        plotPDFfckernewidth_fnc(expSizes,x2_all,x2_b);
    end

    if wtfresamplearefitdist
        % If I use expPD0 to produce a random sample and then I re-estimate new
        % expPD1,2,3,..,nResaples_a then the derived PDFs are accurate? 
        % YES but for this task the best width (~0.5) is different than the best width I
        % found to estimate expPD0 (~3.5)
        figure 
        hold on
        plot(expSizes,expPDF_at_expSizes,'r-x');%,'linewidth',2)
        n = 10000;
        nResaples_a = 1;
        colA = [0, 1, 0];
        colB = [1, 0, 0];
        colors_a = [linspace(colA(1),colB(1),nResaples_a)', linspace(colA(2),colB(2),nResaples_a)', linspace(colA(3),colB(3),nResaples_a)'];
        for aa = 1:nResaples_a
            artifSizes = random(expPD,n,1);
            artifSizes = sort(artifSizes);
            if aa ==1
                expPDF_at_artifSizes = pdf(expPD,artifSizes);
                plot(artifSizes,expPDF_at_artifSizes,'color',colors_a(aa,:))
            end
            widths = [0.5 1 2 3 4 5 6];
            nResaples_b = length(widths);
            colA = [0, 0, 1];
            colB = [1, 0.5, 0];
            colors_b = [linspace(colA(1),colB(1),nResaples_b)', linspace(colA(2),colB(2),nResaples_b)', linspace(colA(3),colB(3),nResaples_b)'];
            vvs = [];
            for bb = 1:length(widths)
                artifPD = fitdist(artifSizes,'Kernel','Kernel','normal','Width',widths(bb));
                artPDF_at_artifSizes = pdf(artifPD,artifSizes);
                vv = plot(artifSizes,artPDF_at_artifSizes,'color',colors_b(bb,:));
                vvs = [vvs vv];
            end
            legend(vvs',char(strcat({'w = '},cellstr(string(widths)))))
        end
    end

    % a] Estimating confidence intervals for the NDF (N_qL, N_qH)
    nboot = 10^4; % number of resamples
    filenameBoot_NDF = strcat('./outputs/data/bootMat_NDF_at_x2_all_fitdist_cutPre', char(string(cutInitialNDFbeforeEstim)),'_kerWidth',num2str(kernelWidth),'_',num2str(nCells),'cells_',sC,'.mat');
    if isfile(filenameBoot_NDF) 
        load(filenameBoot_NDF,'bootM_NDF');
    else
        bootM_NDF = bootstrp(nboot,@(sampledSizes) compute_pd_NDF_PDF(sampledSizes,x2_b,[],kernelWidth,{'NDF'}),expSizes);
        save(filenameBoot_NDF,'bootM_NDF');
    end
    
    % b] Estimating confidence intervals for the PDF (PDF_qL, PDF_qH)
    filenameBoot_PDF = strcat('./outputs/data/bootMat_PDF_at_x_fitdist_cutPre', char(string(cutInitialNDFbeforeEstim)),'_kerWidth',num2str(kernelWidth),'_',num2str(nCells),'cells_',sC,'.mat');
    if isfile(filenameBoot_PDF) 
        load(filenameBoot_PDF,'bootM_PDF');
    else
        xforPDF = x2_all;
        bootM_PDF = bootstrp(nboot,@(sampledSizes) compute_pd_NDF_PDF(sampledSizes,[],xforPDF,kernelWidth,{'PDF'}),expSizes);
        save(filenameBoot_PDF,'bootM_PDF');
    end
    N_qL = arrayfun(@(i) quantile(bootM_NDF(:,i),0.025) , 1:size(bootM_NDF,2), 'UniformOutput',true)';
    N_qH = arrayfun(@(i) quantile(bootM_NDF(:,i),0.975) , 1:size(bootM_NDF,2), 'UniformOutput',true)';
    N_m = mean(bootM_NDF);
    
    PDF_qL = arrayfun(@(i) quantile(bootM_PDF(:,i),0.025) , 1:size(bootM_PDF,2), 'UniformOutput',true)';
    PDF_qH = arrayfun(@(i) quantile(bootM_PDF(:,i),0.975) , 1:size(bootM_PDF,2), 'UniformOutput',true)';
    PDF_m = mean(bootM_PDF);

    save(['./tmpOut_' num2str(cutInitialNDFbeforeEstim) num2str(cutInitialNDFafterEstim) '.mat'])
    
    if cutInitialNDFafterEstim
        [x2_all,x2_b,x2_lb,x2_ub,x2_widths,x2Min,x2Max,N_qL,N_qH,N_data,N_m,xforPDF,PDF_qL,PDF_qH,PDF_data,PDF_m] = cutDistr(x2_all,x2SCell,x2_b,x2_lb,x2_ub,x2_widths,N_qL,N_qH,N_data,N_m,xforPDF,PDF_qL,PDF_qH,PDF_data,PDF_m);
    end
 otherwise
    error('specify how to compute the NDF inside PBEsettings')
end

% Rescale NDF to one
mustintegrateTo = 1; 
% N is the number of particles in the system, i.e. the zeroth moment.
% If I want it to represent a pdf the zeroth moment should be 1->the
% sum of N should be 1
% only N_data is a PDF! The other 2 are resized
% by the same proportionality. If I didn't do this I would have that in
% some points the curves intersect
[N_data,proportionalityK] = rescaleN_st_integralEquals(N_data,mustintegrateTo);
% NB I need to make it so that these two are not PDFs!
N_qL = N_qL*proportionalityK;
N_qH = N_qH*proportionalityK;
N_m = N_m*proportionalityK;

checkDontIntersect(N_qL,N_qH,N_data,PDF_qL,PDF_qH,PDF_data)

%% Compute N_0
% x2bar = sum(N_data.*x2_all);
% dajeoh = P.InitSize_all/x2bar;
% P.N_tot_0_all.(sp) = P.InitSize_all/x2_bar;


%% PLOTS
[minValue,closestIndex] = min(abs(x2_all-max(expSizes)*cutAt_x2Max));
if excludeFirstPoint
 x2_ind_toPlot = 2:closestIndex;
else
 x2_ind_toPlot = 1:closestIndex;
end

% - Plot NDF
figure
plot(x2_all(x2_ind_toPlot),N_data(x2_ind_toPlot),'linewidth',1.5)
hold on
% for nn = 1:100
% plot(x2_all(x2_ind_toPlot),bootM_NDF(nn,x2_ind_toPlot),'k-','linewidth',0.05)
% end
plot(x2_all(x2_ind_toPlot),N_m(x2_ind_toPlot),'--')
plot(x2_all(x2_ind_toPlot),N_qL(x2_ind_toPlot),'blue')
plot(x2_all(x2_ind_toPlot),N_qH(x2_ind_toPlot),'red')
title('estimated N ( = discretized NDF)')

reallyEqual(sum(N_data),1,' = ','Should be the same')

% - Plot PDF
figure
plot(x2_all(x2_ind_toPlot),PDF_data(x2_ind_toPlot),'linewidth',1.5)
hold on
plot(x2_all(x2_ind_toPlot),PDF_m(x2_ind_toPlot),'--')
plot(x2_all(x2_ind_toPlot),PDF_qL(x2_ind_toPlot),'blue')
plot(x2_all(x2_ind_toPlot),PDF_qH(x2_ind_toPlot),'red')
title('estimated PDF')

hold off
reallyEqual(sum(N_data),1,' = ','Should be the same')

% - Plot eCDF
if checkECDF
 eCDF1 = compute_eCDF(expSizes,expSizes);
 eCDF2 = compute_eCDF(expSizes,x2_all);
 plot(expSizes,eCDF1)
 hold on
 plot(x2_all,eCDF2)
 legend({'at exper expSizes','at x2-all'})
 title('eCDF depending on evaluation points')
end

%% define coordinate statistics for the low_biom case
if strcmp(IC_at,'startupScenario')
 aSt.(sp).mu_x2 = aSt.(sp).diamSCell*2;
 stDev_over_mu = 1;% gives me the proportion of stdev/mu
 aSt.(sp).sigma_x2 = aSt.(sp).mu_x2*stDev_over_mu;
end

%% STORE
aSt.(sp).x2_ind_toPlot = x2_ind_toPlot;

aSt.(sp).N_data = N_data(:); 
aSt.(sp).N_qL = N_qL(:); 
aSt.(sp).N_qH = N_qH(:);
aSt.(sp).PDF_data = PDF_data(:); 
aSt.(sp).PDF_qL = PDF_qL(:); 
aSt.(sp).PDF_qH = PDF_qH(:);

aSt.(sp).x2Min = x2Min;
aSt.(sp).x2Max = x2Max;
aSt.(sp).x2SCell = x2SCell;
aSt.(sp).x2_all = x2_all(:);
aSt.(sp).x2_b = x2_b(:);
aSt.(sp).x2_lb = x2_lb(:);
aSt.(sp).x2_ub = x2_ub(:);
aSt.(sp).x2_widths = x2_widths(:);
aSt.(sp) = computeDiam(aSt.(sp));
switch hotToComputeNDF
 case 'ksdensity'
    aSt.(sp).x2_b_subint = x2_b_subint(:);
    aSt.(sp).x2_subint_widths_M = x2_subint_widths_M;
end

%% define weights
aSt = defineObjWeights(aSt,sp);
end
return

function s = areaToVol(a)
 s = ((4/3).*pi.*(sqrt(a/pi)).^3);
return

function d_all = diam(v_all)
 d_all = 2*((3/(4*pi))*v_all).^(1/3); 
return

function aSt = computeDiam(aSt)
 aSt.d_all = diam(aSt.x2_all);
 aSt.d_lb = diam(aSt.x2_lb);
 aSt.d_ub = diam(aSt.x2_ub);
return

function [x2_all,x2_b,x2_lb,x2_ub,x2_widths,x2Min,x2Max,varargout] = cutDistr(x2_all,x2SCell,x2_b,x2_lb,x2_ub,x2_widths,varargin)
 % CUT. Once the NDF has been obtained, I cut out the distr < x2SCell
 global nCells
 x2_all = x2_all(end-nCells+1:end);
 if ~reallyEqual(x2_all(1),x2SCell)
    error('My a-posteriori cut is not working properly')
 end
 x2_b = x2_b(end-nCells:end);
 x2_lb = x2_lb(end-nCells+1:end);
 x2_ub = x2_ub(end-nCells+1:end);
 
 % IMPO min and max of the coordinate space are boundary values not
 % grid points, otherwise I have problems of negative numbers when I
 % normalize between x2Min = 0 - 1=x2Max
 x2Min = min(x2_b);
 x2Max = max(x2_b);
 x2_widths = x2_widths(end-nCells+1:end);
 if nargin>= 10
    N_qL = varargin{1};
    N_qH = varargin{2};
    N_data = varargin{3};
    N_m = varargin{4};
    N_qL = N_qL(end-nCells+1:end);
    N_qH = N_qH(end-nCells+1:end);
    N_data = N_data(end-nCells+1:end); 
    N_m = N_m(end-nCells+1:end); 
    varargout{1} = N_qL;
    varargout{2} = N_qH;
    varargout{3} = N_data;
    varargout{4} = N_m;
 end
 if nargin ==15
    xforPDF = varargin{5};
    idxToKeep = find(xforPDF>x2Min & xforPDF<x2Max);
    xforPDF = xforPDF(idxToKeep);
    PDF_qL = varargin{6};
    PDF_qH = varargin{7};
    PDF_data = varargin{8};
    PDF_m = varargin{9};
    PDF_qL = PDF_qL(idxToKeep);
    PDF_qH = PDF_qH(idxToKeep);
    PDF_data = PDF_data(idxToKeep); 
    PDF_m = PDF_m(idxToKeep); 
    varargout{5} = xforPDF;
    varargout{6} = PDF_qL;
    varargout{7} = PDF_qH;
    varargout{8} = PDF_data;
    varargout{9} = PDF_m;
 end
return

function aSt = defineObjWeights(aSt,sp)
% - Defining the weights that weigh the difference data-simResult in the
% obj function: (yOneParam-yData).*objWeights. Defined inversely
% proportional to the measure of uncertainty (widht of the confidence interval)
 aSt.(sp).weights.equal = ones(size(aSt.(sp).N_qH));
 aSt.(sp).weights.NDF = 1./(aSt.(sp).N_qH-aSt.(sp).N_qL);
 aSt.(sp).weights.PDF = 1./(aSt.(sp).PDF_qH-aSt.(sp).PDF_qL);
return

function [] = checkDontIntersect(N_qL,N_qH,N_data,PDF_qL,PDF_qH,PDF_data)
if ~all(N_qL<= N_data) || ~all(N_data<=N_qH)
 idxL = find(N_qL>N_data);
 tuttOK = true;
 for yy = 1:length(idxL)
    if ~reallyEqual(N_qL(idxL(yy)),N_data(idxL(yy)))
        tuttOK = false;
        break;
    end
 end
 idxH = find(N_qH<N_data);
 for yy = 1:length(idxH)
    if ~reallyEqual(N_qH(idxH(yy)),N_data(idxH(yy)))
        tuttOK = false;
        break;
    end
 end
 if ~tuttOK
    error('the lower and upper confidence intervals should be on the 2 sides of the average. It is possible that I need to reduce the kernelWidth');
 end
 error('the lower and upper confidence intervals should be on the 2 sides of the average. It is possible that I need to reduce the kernelWidth');
end
if ~all(PDF_qL<= PDF_data) || ~all(PDF_data<=PDF_qH)
 idxL = find(PDF_qL>PDF_data);
 tuttOK = true;
 for yy = 1:length(idxL)
    if ~reallyEqual(PDF_data(idxL(yy)),PDF_qL(idxL(yy)))
        tuttOK = false;
        break;
    end
 end
 idxH = find(PDF_qH<PDF_data);
 for yy = 1:length(idxH)
    if ~reallyEqual(PDF_data(idxH(yy)),PDF_qH(idxH(yy)))
        tuttOK = false;
        break;
    end
 end
 if ~tuttOK
    error('the lower and upper confidence intervals should be on the 2 sides of the average. It is possible that I need to reduce the kernelWidth');
 end
end
return





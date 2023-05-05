function [] = checkNDFtoPDF(expSizes,N_res,x2_all,x2_b,x2_lb,x2_ub,x2_widths)
%% AD and CVM are based on CDF -> Check that the following shapes are close:
% a] expSizes -*> expPD -> expCDF -> NDF -> mysampling -*> compPD -> compPDF 
% b] expSizes -> expPD -> expPDF
% Conclusion: are close but also separate = > for the experimental
% distribution I consider it to be equal to the first tPoint of compPDF,
% instead of expPDF

global kernelWidth

% IMPO: here I am only considering all the coordinate values x2_ind_toPlot = 1:size(N_res,2);
oneTpoint = 1;
oneNDF = N_res(oneTpoint,:);

% prendo N_res -> create its "sample" vector S: if at x1 NDF = 2.33456 and
% nDecimalsToInclude = 2 → I put 233 times the element x1 in S. I do the
% same for all grid points 
nDecimalsToInclude = 4;
nPartStart_all = round(oneNDF*10^nDecimalsToInclude);
S = [];
for i = 1:length(oneNDF)
 % take the number of particles to be sampled for each grid interval/cell
 nPartStartThisInt = nPartStart_all(i);
 if nPartStartThisInt>0
    % random sampling n-2 particles. The last 2 are added to make
    % the meanThisSample on point
    unifSampleThisInt = x2_lb(i) + x2_widths(i)*rand(1,nPartStartThisInt-2);
    sampleThisInt = unifSampleThisInt;
    nThisSample = length(sampleThisInt);
    meanThisSample = mean(sampleThisInt);
    gridPoint = x2_all(i);
    % I add n new particles as long as the grid point is not
    % representative of the whole interval
    nP = 1;
    while nThisSample<nPartStartThisInt 
        fprintf('n = %g particle added to int %g \n',nP,i)
        % adjust the size of the new particle so that the mean in the
        % interval matches
        newPart = gridPoint*(nThisSample+1)-sum(sampleThisInt);
        % If newPart falls out of the boundaries I make it so that it lays
        % inbetween (not exactly on the boundary, just in case)
        if newPart<x2_lb(i)
            newPart = x2_lb(i)+x2_widths(i)/100;
        end
        if newPart>x2_ub(i)
            newPart = x2_ub(i)-x2_widths(i)/100;
        end
        sampleThisInt = [sampleThisInt newPart];
        nThisSample = length(sampleThisInt);
        meanThisSample = mean(sampleThisInt);
        nP = nP+1;
    end
    % add sampleThisInt to the main sample
    S = [S sampleThisInt];
 end
end
S = sort(S(:));

% Check that the following shapes are close:
% a] expSizes -*> expPD -> expCDF -> NDF -> mysampling -*> compPD -> compPDF 
% b] expSizes -> expPD -> expPDF
% I would expect similar widths to be needed in the 2 fitdist processes '-*>'
figure
% plot(x2_all,onePDF/max(onePDF),'b*-','linewidth',2)
% fitdist → pd object that resembles the shape of the NDF
widths = sort([0.1:0.2:3 kernelWidth]);
red = [1, 0, 0];
green = [0, 1, 0];
colors_p = [linspace(red(1),green(1),length(widths))', linspace(red(2),green(2),length(widths))', linspace(red(3),green(3),length(widths))'];
hold on
hhs = [];
for ii = 1:length(widths)
 aWi = widths(ii);
 compPDbase = fitdist(S,'Kernel','Kernel','normal','Width',aWi);
 compPDF = pdf(compPDbase,x2_all);
 if aWi ==kernelWidth
    plot(x2_all,compPDF,':','color',colors_p(ii,:),'linewidth',1.5)
    compPDF = pdf(compPDbase,expSizes);
    hh = plot(expSizes,compPDF,'color',colors_p(ii,:),'linewidth',1.5);
 else
    plot(x2_all,compPDF,':','color',colors_p(ii,:));
    compPDF = pdf(compPDbase,expSizes);
    hh = plot(expSizes,compPDF,'color',colors_p(ii,:));
 end
 hhs = [hhs hh];
end 
expPD = compute_pd_NDF_PDF(expSizes,[],[],kernelWidth,{'pd'});
compPDF = pdf(expPD,expSizes);
hh = plot(expSizes,compPDF,'color','b','linewidth',1.5);
hhs = [hhs hh];
legend(hhs,[ strcat('compPD w = ',cellstr(string(widths))) {'expPDF'}]);
hold off


% dist object expPD coming from the experimental sample -> what if I
% use it to evaluate the pdf at different locations than the ones used
% for the fitting? I would expect them to be on top of the line
% make sure x values do not fall outside of the range
expPD = fitdist(expSizes,'Kernel','Kernel','normal','Width',kernelWidth);
expPDF_at_expSizes = pdf(expPD,expSizes);
hold on
plot(expSizes,expPDF_at_expSizes,'r-x');%,'linewidth',2)
expPDF_at_x2_all = pdf(expPD,x2_all);
expPDF_at_x2_b = pdf(expPD,x2_b);
plot(x2_all,expPDF_at_x2_all,'g-o');%,'linewidth',2)
plot(x2_b,expPDF_at_x2_b,'b-*');%,'linewidth',2)
% plotPDFfckernewidth_fnc(expSizes,x2_all,x2_b)
return

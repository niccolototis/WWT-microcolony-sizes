function [aSt,P] = PBEParams(P,plotsVisible)
global initialNDF 	InitBiom_all  	 TRUE RESC0to1  normalizeCoordinates  
global nCells sC mu_adh_ sigma_adh_ 	 

% - If I am using different parameters from 'baseling ' control that ASM1seq
% parameters are the same as the one used here in the PBE
if strcmp(P.parSetup,'chosenSetASM1seq')
 nVars = length(P.parToChange);
 for i = 1:nVars
    reallyEqual(P.(char(P.parToChange(i))),P.NewValue(i),' = ','The parameter values used in the ASM1seq model are overwritten somewhere');
 end
end

% -- Biologically relevant parameters
mu_adh_ = [5,5,5];
sigma_adh_ = [0.01,0.01,0.01];
P.InitBiom_all = InitBiom_all; % defined as a global variable inside defineInitialConditions_t0 so that I dont need to return P as an output argument
P.InitSize_all = P.InitBiom_all./P.sludgeDensity; % this is in m^3 
% P.InitSize_all = P.InitSize_all*1e18; % this is in micron^3 

if strcmp(initialNDF,'x2MarginalByKDEstimate')
 TRUE = define_mu_sigma_x1(TRUE,P);
 [TRUE,P] = kdestimation(TRUE,P);
else
 TRUE = [];
 TRUE = define_mu_sigma_x1(TRUE,P);
 TRUE = define_mu_sigma_x2(TRUE,P);
 TRUE = defineGrid_x2(TRUE,sC,nCells);
end
if plotsVisible
 plotCoordinates(TRUE)
end
aSt = TRUE; RESC0to1=[]; 
if normalizeCoordinates
 RESC0to1 = normalizeCoord_x1(TRUE,RESC0to1);
 RESC0to1 = normalizeCoord_x2(TRUE,RESC0to1);
 if strcmp(initialNDF,'x2MarginalByKDEstimate')
    RESC0to1 = normalizeGrid_x2(TRUE,RESC0to1);
 else
    RESC0to1 = defineGrid_x2(RESC0to1,sC,nCells);
 end
 if plotsVisible
    plotCoordinates(RESC0to1)
 end
 aSt = RESC0to1;
end
aSt = define_PDF_PVFs(aSt,plotsVisible);
aSt = addForGrowth(aSt);
return

function TRUE = define_mu_sigma_x1(TRUE,P)
global sigma_adh_ mu_adh_ 
% Define TRUE values (not normalized)
for aSpPBE = P.speciesPBE
sp = char(P.spNames{aSpPBE}); 
 % 1) for x1: mu_adh, dev_st_adh, sigma_adh
 TRUE.(sp).mu_x1 = mu_adh_(aSpPBE);
 TRUE.(sp).sigma_x1 = sigma_adh_(aSpPBE);
 TRUE.(sp).x1Max = TRUE.(sp).mu_x1*4; %@param
 TRUE.(sp).x1Min = 0;
end
return

function TRUE = define_mu_sigma_x2(TRUE,P)
global nCells
% Define TRUE values (not normalized)
for aSpPBE = P.speciesPBE
sp = char(P.spNames{aSpPBE}); 
 % 2) for x2: mu_v, dev_st_v, sigma_v (+ [vMin vMax], InitVol_whole_aSp, InitBiom_aSp, N_tot(initial number of particles))% Particle volume coordinate
 TRUE.(sp).diamSCell = 1; % [micron] diameter single cell. This is the diameter of the minimal particle that can be eroded --> first grid point is given this value
 TRUE.(sp).diamAveragePart = 13; % [micron] this is the diameter of the average microcolony
 TRUE.(sp).diamMaxMiCol = 100; % [micron] ; 150e-06 m
 TRUE.(sp).mu_x2 = diamToVol(TRUE.(sp).diamAveragePart); % In this case they are all equal
 TRUE.(sp).x2SCell = diamToVol(TRUE.(sp).diamSCell);
 TRUE.(sp).x2Max = diamToVol(TRUE.(sp).diamMaxMiCol); 
 stDev_over_mu = 1/5;% gives me the proportion of stdev/mu
 TRUE.(sp).sigma_x2 = TRUE.(sp).mu_x2*stDev_over_mu;
 %-- create mu 2x1 and sigma 2x2 matrices later used to build the
 %probability density function in defineIC
 TRUE.(sp).mu = [TRUE.(sp).mu_x1 TRUE.(sp).mu_x2]; % in cm
 TRUE.(sp).sigma = [TRUE.(sp).sigma_x1 TRUE.(sp).sigma_x2].*eye(2); % in cm
 
 TRUE.(sp).x2_ind_toPlot = 1:nCells;
end
return

function aSt = defineGrid_x2(aSt,sC,nCells) 
global speciesPBE spNames logBasis 
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE});
aStp = aSt.(sp);
switch sC
 case 'lin' 
    % x2MmiCol and x2SCell are set at the grid points of the first and
    % last cell respectively. x2Min is set to the lower bound of the
    % first cell
    % --
    % I first create (nCells-1) cells between x2SCell and x2MmiCol
    % then I shift them to the left to have x2SCell centered at the
    % first grid point
    % then I add the last cell at the right end side
    % NB: SCell refers to a single bacterial cell. cellWidth refers to
    % the cell of the discretization grid (interval among 2 boundaries)
    aStp.x2_b = linspace(aStp.x2SCell,aStp.x2Max,nCells)';
    cellWidth = (aStp.x2Max-aStp.x2SCell)/(nCells-1);
    aStp.x2_b = aStp.x2_b-cellWidth/2;
    aStp.x2_b = [aStp.x2_b ; aStp.x2_b(end)+cellWidth];
    aStp.x2Min = aStp.x2_b(1);   
    aStp.x2Max = aStp.x2_b(end); 
    x2_lb = aStp.x2_b(1:end-1);
    x2_ub = aStp.x2_b(2:end);
    aStp.x2_all = (x2_lb+x2_ub)/2; % reference at the midpoint
    MEstr = 'I want the first grid point to be centred at the single cell volume and the last to be centered at the MaxMiCol volume';
    [~,A,~] = reallyEqual(aStp.x2_all(1),aStp.x2SCell,'~',MEstr);
    aStp.x2_all(1) = A;
    aStp.x2SCell = A;
    
    aStp.x2_b(1) = aStp.x2SCell;
    aStp.x2Min = aStp.x2_b(1);   
    aStp.x2_lb = aStp.x2_b(1:end-1);
    aStp.x2_ub = aStp.x2_b(2:end);
    
    [~,A,~] = reallyEqual(aStp.x2_all(end),aStp.x2Max,'~',MEstr);
    aStp.x2_all(end) = A;
    aStp.x2Max = A;
 case 'log'       
    expSCell = log(aStp.x2SCell)/log(logBasis);
    expMaxMCol = log(aStp.x2Max)/log(logBasis);
    exp_b = linspace(expSCell,expMaxMCol,nCells)';
    expCellWidth = (expMaxMCol-expSCell)/(nCells-1);
    exp_b = exp_b-expCellWidth/2;
    exp_b = [exp_b ; exp_b(end)+expCellWidth];
    exp_lb = exp_b(1:end-1);
    exp_ub = exp_b(2:end);
    exp_all = (exp_lb+exp_ub)/2;
    if ~strcmp(string(exp_all(1)),string(expSCell)) || ~strcmp(string(exp_all(end)),string(expMaxMCol)) 
        ME = MException('I want the exp of the first grid point to be centred at the exp of single cell (expSCell) and the exp of the last cell to be centered at the expMaxMCol');
        throw(ME);
    end
    aStp.x2_b = logBasis.^exp_b;
    aStp.x2_lb = aStp.x2_b(1:end-1);
    aStp.x2_ub = aStp.x2_b(2:end);
    aStp.x2Min = aStp.x2_b(1);   
    aStp.x2Max = aStp.x2_b(end); 
    aStp.x2_all = logBasis.^exp_all;
end
aStp = computeDiam(aStp);
aSt.(sp) = aStp;
end
return

function RESC0to1 = normalizeCoord_x1(TRUE,RESC0to1)
global spNames speciesPBE
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE});
 RESC0to1.(sp).mu_x1 = TRUE.(sp).mu_x1/TRUE.(sp).x1Max;
 RESC0to1.(sp).sigma_x1 = TRUE.(sp).sigma_x1/TRUE.(sp).x1Max;
 RESC0to1.(sp).x1Min = TRUE.(sp).x1Min/TRUE.(sp).x1Max;
 RESC0to1.(sp).x1Max = TRUE.(sp).x1Max/TRUE.(sp).x1Max; %@NB this has to result=1!   
end
return

function RESC0to1 = normalizeCoord_x2(TRUE,RESC0to1)
global spNames speciesPBE
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE}); 
 scalingFac = TRUE.(sp).x2Max;
 RESC0to1.(sp).mu_x2 = TRUE.(sp).mu_x2/scalingFac;
 RESC0to1.(sp).sigma_x2 = TRUE.(sp).sigma_x2/scalingFac;
 RESC0to1.(sp).x2SCell = TRUE.(sp).x2SCell/scalingFac;
 RESC0to1.(sp).x2Max = TRUE.(sp).x2Max/scalingFac; %@NB this has to result=1!  
 RESC0to1.(sp).x2Min = TRUE.(sp).x2Min/scalingFac;
 RESC0to1.(sp).x2Max = TRUE.(sp).x2Max/scalingFac; %@NB this has to result=1!  
 %-- create mu 2x1 and sigma 2x2 matrices later used to build the
 %probability density function in defineIC
 RESC0to1.(sp).mu = [RESC0to1.(sp).mu_x1 RESC0to1.(sp).mu_x2]; % in cm
 RESC0to1.(sp).sigma = [RESC0to1.(sp).sigma_x1 RESC0to1.(sp).sigma_x2].*eye(2); % in cm
end
return

function RESC0to1 = normalizeGrid_x2(TRUE,RESC0to1)
global spNames speciesPBE initialNDF hotToComputeNDF 
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE}); 
 scalingFac = TRUE.(sp).x2Max;
 RESC0to1.(sp).x2_all_TRUE_micronCube = TRUE.(sp).x2_all; % Stored for later plotting. Is in cubic micron
 RESC0to1.(sp).x2_b_TRUE_micronCube = TRUE.(sp).x2_b;
 RESC0to1.(sp).x2SCell_TRUE_micronCube = TRUE.(sp).x2SCell;
 RESC0to1.(sp).x2_widths_TRUE_micronCube = TRUE.(sp).x2_widths; % Stored for later plotting. Is in cubic micron
 RESC0to1.(sp).x2_all = TRUE.(sp).x2_all/scalingFac;
 RESC0to1.(sp).x2_lb = TRUE.(sp).x2_lb/scalingFac;
 RESC0to1.(sp).x2_ub = TRUE.(sp).x2_ub/scalingFac;
 RESC0to1.(sp).x2_b = TRUE.(sp).x2_b/scalingFac;
 RESC0to1.(sp).x2_widths = TRUE.(sp).x2_widths/scalingFac;
 RESC0to1.(sp).x2SCell = TRUE.(sp).x2SCell/scalingFac;
 RESC0to1.(sp).x2Max = TRUE.(sp).x2Max/scalingFac; %@NB this has to result=1!  
 RESC0to1.(sp).x2Min = TRUE.(sp).x2Min/scalingFac;
 RESC0to1.(sp).x2Max = TRUE.(sp).x2Max/scalingFac; %@NB this has to result=1!  
 if strcmp(initialNDF,'x2MarginalByKDEstimate') 
    if strcmp(hotToComputeNDF,'x2MarginalByKDEstimate') 
        RESC0to1.(sp).x2_b_subint = TRUE.(sp).x2_b_subint/scalingFac;
        RESC0to1.(sp).x2_subint_widths_M = TRUE.(sp).x2_subint_widths_M/scalingFac;
    end
    % - As I change the x coordinates, the pdf function (the y
    % coordinate) changes, but the zeroth moments do not! So I do not
    % need to rescale them       
    RESC0to1.(sp).N_data = TRUE.(sp).N_data;
    RESC0to1.(sp).N_qL = TRUE.(sp).N_qL;
    RESC0to1.(sp).N_qH = TRUE.(sp).N_qH;
    
    % - indexes
    RESC0to1.(sp).x2_ind_toPlot = TRUE.(sp).x2_ind_toPlot;
 end
end
return

function [] = plotCoordinates(aSt)
global spNames speciesPBE initialNDF 
figure
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE}); 
 % 1) for x1: mu_adh, dev_st_adh, sigma_adh
 pd = makedist('Normal','mu',aSt.(sp).mu_x1,'sigma',aSt.(sp).sigma_x1);
 x = linspace(max(0,(aSt.(sp).mu_x1-3*aSt.(sp).sigma_x1)),(aSt.(sp).mu_x1+3*aSt.(sp).sigma_x1),50);
 y = pdf(pd,x);
 if ~any(y) % if y is all zeros
    ME = MException('It is probably due to the fact that sigma is too small');
    throw(ME)
 end
 plot(x,y)
 hold on
end
hold off
title('Adh distribution')
legend(spNames(speciesPBE));

figure
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE}); 
 if strcmp(initialNDF,'x2MarginalByKDEstimate') 
    x = aSt.(sp).x2_all;
    plot(x,aSt.(sp).N_data)
    hold on
    plot(x,aSt.(sp).N_qL)
    plot(x,aSt.(sp).N_qH)
 else
    x = linspace(max(0,(aSt.(sp).mu_x2-3*aSt.(sp).sigma_x2)),(aSt.(sp).mu_x2+3*aSt.(sp).sigma_x2),50);
    y = normpdf(x,aSt.(sp).mu_x2,aSt.(sp).sigma_x2);
    plot(x,y)
    hold on
 end
 if ~any(y) % if y is all zeros
    ME = MException('It is probably due to the fact that sigma is too small');
    throw(ME)
 end
 plot(aSt.(sp).x2_lb,zeros(size(aSt.(sp).x2_lb)),'bx')
 hold on
 plot(aSt.(sp).x2_ub,zeros(size(aSt.(sp).x2_ub)),'bx')
end
hold off
title('Size distribution')
legend(spNames(speciesPBE));
return


function aStp = computeDiam(aStp)
 aStp.d_all = diam(aStp.x2_all);
 aStp.d_lb = diam(aStp.x2_lb);
 aStp.d_ub = diam(aStp.x2_ub);
return

function d_all = diam(v_all)
 d_all = 2*((3/(4*pi))*v_all).^(1/3); 
return

function vol = diamToVol(diam)
vol = 4/3*pi*(diam/2)^3;
return

function aSt = define_PDF_PVFs(aSt,plotsVisible)
global spNames speciesPBE
for aSpPBE = speciesPBE
 sp = char(spNames{aSpPBE});
 Ss = aSt.(sp);
 % - Compute PDF 
 n_pdf_data = Ss.N_data./Ss.x2_widths_TRUE_micronCube;
 n_pdf_qL = Ss.N_qL./Ss.x2_widths_TRUE_micronCube;
 n_pdf_qH = Ss.N_qH./Ss.x2_widths_TRUE_micronCube;
 % - Compute PVF 
 n_pvf_data = n_pdf_data.*Ss.x2_all_TRUE_micronCube;
 n_pvf_qL = n_pdf_qL.*Ss.x2_all_TRUE_micronCube;
 n_pvf_qH = n_pdf_qH.*Ss.x2_all_TRUE_micronCube;
 
 % - Store in output structure
 aSt.(sp).n_pdf_data = n_pdf_data(:); 
 aSt.(sp).n_pdf_qL = n_pdf_qL(:); 
 aSt.(sp).n_pdf_qH = n_pdf_qH(:); 
 % -
 aSt.(sp).n_pvf_data = n_pvf_data(:); 
 aSt.(sp).n_pvf_qL = n_pvf_qL(:); 
 aSt.(sp).n_pvf_qH = n_pvf_qH(:); 
 
 % - Plot pdf
 x2_ind_toPlot = Ss.x2_ind_toPlot;
 x2_all_TRUE_micronCube = Ss.x2_all_TRUE_micronCube;
 if plotsVisible
    figure
    plot(x2_all_TRUE_micronCube(x2_ind_toPlot),n_pdf_data(x2_ind_toPlot))
    hold on
    plot(x2_all_TRUE_micronCube(x2_ind_toPlot),n_pdf_qL(x2_ind_toPlot),'blue')
    plot(x2_all_TRUE_micronCube(x2_ind_toPlot),n_pdf_qH(x2_ind_toPlot),'red')
    hold off
    title('pdf TRUE')
 end
 
 % - Plot pvf
 x2_ind_toPlot = Ss.x2_ind_toPlot;
 x2_all_TRUE_micronCube = Ss.x2_all;
 if plotsVisible
    figure
    plot(x2_all_TRUE_micronCube(x2_ind_toPlot),n_pvf_data(x2_ind_toPlot))
    hold on
    plot(x2_all_TRUE_micronCube(x2_ind_toPlot),n_pvf_qL(x2_ind_toPlot),'blue')
    plot(x2_all_TRUE_micronCube(x2_ind_toPlot),n_pvf_qH(x2_ind_toPlot),'red')
    hold off
    title('pvf TRUE')
 end
end
return

function aSt = addForGrowth(aSt)
global m_0 spNames speciesPBE
% -- "diffuse particles" for growth
m_0 = 1; % complexive volume/size of of diffuse particles x0
for aSpPBE = speciesPBE
sp = char(spNames{aSpPBE});
 aSt.(sp).x2_x0 = min(aSt.(sp).x2_all(aSt.(sp).x2_all>0))/100; %size 
 aSt.(sp).N_x0 = m_0/aSt.(sp).x2_x0; %number 
end
return








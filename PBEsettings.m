function P = PBEsettings(P)
global thr minNeg kernelWidth cutInitialNDFbeforeEstim   x2MinIsSCell  cutAt_x2Max T_sim plotFrames plotPVF isRunning tspanCut excludeFirstPoint 
global epsylon breakInTwo particlesPileUp  initialNDF IC_at  hotToComputeNDF cutInitialNDFafterEstim  nTrapezPerCell nObj saveOptimalRun
global importanceSampling lostFirstMoment_lb lostFirstMoment_ub tFlaglb tFlagub subsAsInput ModelBiomass out_vAV 
global myPlotOff impulseDistribution HexColors normalizeCoordinates epsyBreak epsyGrow nInconsGrowth nInconsBreak sweep_lambda   
global spNames lastSub  modesGrHET  approx_lb_max approx_ub_max FirstRefNoBreak lambda redistrFn debugParallel lambdasToSweep identify_lambda
global solubPlotted solubNames scaleDots normalized_N_TI_t0 plotsWhileStatTests wtfresamplearefitdist sC logBasis nChar nCells nu adhesins adhDistribution BreakageRate wtbwidthftkerneloexpPD

% ### 1 ###
% - N_PDF_PVF conversions:
%     - kdestimation.m:
%         - ksdensity → produce kerned density estimate (is a pdf)
%         - pdf → N_data that sums to 1
%     - PBEparams.m: N_data → pdf,pvf
%     - DefineIC: N_0 = N_data
%     - DefineObjWeights based on pdf,pvf
%     - runPBEforOneParam: N_res → pdf,pvf
%     - computeObjFcn: pdf,pvf used in the computation
%     - ArrangeToPlot: only the datapoints useful to be plotted are selected

% ### 2 ###
% x2_indToPlot specifies the positions of the distribution I am interested
% to use for plotting
% Thus runPBEforOneParam returns the whole distribution, which is then cut
% with x2_indToPlot both in function prepareToPlot.m (redundant)
 
isRunning = 'PBE';

%% Options for parallel runs
debugParallel = true;

%% Sweep_lambda options
identify_lambda = false; 
% lambdasToSweep = [0.0001 0.001 0.01 0.1 0.5 0.8 1 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10 11 12 13 14 15 23 30 40 50 100];

% aSt.(sp).objVal.all = [...
% aSt.(sp).objVal.notWeighted.N,...
% aSt.(sp).objVal.notWeighted.pdf,...
% aSt.(sp).objVal.notWeighted.pvf,...
% aSt.(sp).objVal.weighted.N,...
% aSt.(sp).objVal.weighted.pdf,...
% aSt.(sp).objVal.weighted.pvf,...
% ];
saveOptimalRun = false;
nObj = 6; % 6 different objective functions evaluations


%% PBE numerical solution scheme options
% NB these two options are different! I can have normalizedNDF_t0 = true,
% normalizedPBE = false to avoid numerical problems for having the NDF too
% large. In this case, the original NDF dynamics are recovered a posteriori
nu = 2;
% -- Parameters of numerical schemes CAT and MOC
nCells = 100; % Discretization points for CAT
switch IC_at 
 case 'SS_computed'
    cutAt_x2Max = 0.05; 
 case 'startupScenario'
    cutAt_x2Max = 0.10; 
end
excludeFirstPoint = false;
nTrapezPerCell = 3; % Used to "integrate" ksdensity => N
% -
sC = 'log';% Scale of discretization grid %'log';%'lin'
logBasis = 2;
% - MOC
nChar = 1; % Number of characteristic curves for MOC

%% Initial distribution options - check kdestimation.m
cutInitialNDFbeforeEstim = true;
cutInitialNDFafterEstim = false;
if (cutInitialNDFbeforeEstim && cutInitialNDFafterEstim) || (~cutInitialNDFbeforeEstim && ~cutInitialNDFafterEstim)
 error('cannot be both true or both false')
end
hotToComputeNDF = 'fitdist'; % {'fitdist','ksdensity'} 
kernelWidth = 1;
wtbwidthftkerneloexpPD = false; % whats the best width for the kernel of the expPD
wtfresamplearefitdist = false; % if I use expPD0 to produce a random sample and then I re-estimate new
normalized_N_TI_t0 = true; 
initialNDF = 'x2MarginalByKDEstimate';%'normal_2D';%'impulse_2D';
FirstRefNoBreak = true; % Once I have checked that the first reference is at x2SCell I can activate this, later used in intBrKerWWT and odeFunPBE
x2MinIsSCell = true;
adhDistribution = 'normal'; %'lognormal';%normal'; %'lognormal';

%% Options for when I compute my similarity tests
plotsWhileStatTests = false;

%% Simulation options
adhesins = 1;
BreakageRate = 'linear';
lambda = 10; % Breakage 
redistrFn = 'U-shape'; %'constant'; %

breakInTwo = true;
importanceSampling = false;
impulseDistribution = false;
ModelBiomass = false; %%%%%%%% DOUBLE CHECK %If true the x2 coordinate is expressed as Kg and not as volume (m^3)
normalizeCoordinates = true;
myPlotOff = false;
particlesPileUp = false; approx_lb_max=0; approx_ub_max=0;
% control pqrqmeters
[epsyBreak, epsyGrow, nInconsGrowth, nInconsBreak,lostFirstMoment_lb, lostFirstMoment_ub, tFlaglb, tFlagub,out_vAV] = deal(0);
thr = 0; % when variables go below this thr, they are reset to thr within the odefun
epsylon = 0; % defines the allowed imprecision when total births - total deaths should be the same. NUmerical errors. Is this numerical diffusion?
minNeg = 0;


%% Time simulation parameters
% - Here I consider a shorter tspan just to lower computational effort.
% After a certain time I saw that the PSD does not change. 
% To try to reduce the computational time I can reduce the tspan of the
% actual simulation by setting P.T_sim_PBE. The most accurate tspan is 50 days. After the first movies
% however I see that from day 30-35 to day 50 the PSD is already pretty
% much stationary. So I set P.T_sim_PBE to 30 or 35, then I take the final
% PSD as the steady state PSD and I sort of create artificial data from 35
% to 50 by coy-paste the same PSD
tspanCut = false;
T_sim_PBE = 30; %days. It should be a multiple of 0.5 = operation cycle length 

if T_sim<= T_sim_PBE
 tspanCut = false;
end
if P.tspanCut
 lastFragNo = [];
 for nf = 1:length(fieldnames(P.tspans_ASM1seq))
    oneFrag = P.tspans_ASM1seq.(strcat(['f',num2str(nf)]));
    if reallyEqual(P.T_sim_PBE,oneFrag(end))
        lastFragNo = nf;
        P.nPhases_all_ASM1seq = P.nPhases_all;
        P.nPhases_all = nf;
        break;
    end
 end
 if isempty(lastFragNo)
    ME = MException('Something going wrong I should have found a fragment which ends at this time');
    throw(ME);
 end
end

%% Plotting options
% Plotting the Probability Volume Function PVF, like probability mass
% function. Is N(i.e. zeroth moment) in each cell of the grid times the
% volume(size) of the reference point
plotPVF = true;

% Reduce the number of time points that are retained in matrices after PBE
% simulation
plotFrames = 100;
if P.ntpts_ASM1seq>1000 
 Q = floor(P.ntpts_ASM1seq/plotFrames);
 P.idx_tpts_plot_PBE = 1:Q:P.ntpts_ASM1seq;
 P.idx_tpts_plot_PBE = [P.idx_tpts_plot_PBE P.ntpts_ASM1seq];
 P.ntpts_plot_PBE = length(P.idx_tpts_plot_PBE);
else
 P.ntpts_plot_PBE = P.ntpts_ASM1seq;
 P.idx_tpts_plot_PBE = 1:P.ntpts_plot_PBE;
end

% HET = 1;AOB=2;NOB=3;
HexColors = {'00ba38','619cff','ca7f28'};

C = linspecer(P.nSp); 
ll = 1;
for aSp = P.speciesModeled
 sp = char(spNames{aSp});
 spColors.(sp) = C(ll,:);
 ll = ll+1;
end
scaleDots = 'bo';
solubPlotted = P.solubles;
lastSub = solubPlotted(end);
solubNames = P.solNames(solubPlotted);
modesGrHET = {'AER', 'NO2','NO3'};

%% Put all global variables into structure

namesGlobal = who('global'); %# A cell array of variable names
namesLocal = who(); 
namesGLocal = intersect(namesGlobal,namesLocal); 
for iVar = 1:numel(namesGLocal)
 P.(namesGLocal{iVar}) = eval(namesGLocal{iVar}); % [EDITED]
end

return

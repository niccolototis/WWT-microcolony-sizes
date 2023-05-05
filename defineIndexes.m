function P = defineIndexes(P)
global sI sO2 sNO sND sNH4 sALK sS xI xS xP xND nPart nSolub solubles particles speciesPBE speciesUnif nSpPBE nSpU
global solNames solLines partLines partNames spNames nSp HET AUT speciesModeled folderInput folderOutput folderPlots T_sim rqrmtsIdxs
global tCOD tN MLVSS OUR rqrmNames IC_at biom_lower_than_SS computeSS_for_later_IC biomFract0
 HET = 1; AUT=2; xI=3; xS=4; xP=5; xND=6;
 sI = 1; sS=2; sO2=3; sNH4=4; sNO=5; sND=6; %sALK=7; 
 tCOD = 1;tN=2;MLVSS=3;OUR=4;
 spNames = {'HET','AUT'};
 nSp = length(spNames);
 particles = [HET AUT xI xS xP xND];
% solubles = [sI sS sO2 sNH4 sNO sND sALK];
 solubles = [sI sS sO2 sNH4 sNO sND];
 rqrmtsIdxs = [tCOD,tN,MLVSS,OUR];

 nPart = length(particles);
 nSolub = length(solubles);
% solNames = {'I','S','O2','NH4','NO','ND','ALK'};
 solNames = {'I','S','O2','NH4','NO','ND'};
 partNames = {'xI','xS','xP','xND'};
 rqrmNames = {'tCOD','tN','MLVSS','OUR'};
 [solLines{1:length(solNames)}, partLines{1:length(partNames)}] = deal('-');

 %-
 speciesPBE = AUT; %[AOB NOB]; % for these species the PSD is taken into account
 speciesUnif = HET; % These species are represented with just a uniform biomass variable, no size attribute
 speciesModeled = [speciesUnif speciesPBE];
 nSpPBE = length(speciesPBE);
 nSpU = length(speciesUnif);
 
 
if ~exist(['./outputs/data/'], 'dir')
   mkdir(['./outputs/data/'])
end
if ~exist(['./outputs/parameter_modifyOnce/'], 'dir')
   mkdir(['./outputs/parameter_modifyOnce/'])
end
if ~exist(['./outputs/figures/'], 'dir')
   mkdir(['./outputs/figures/'])
end
if ~exist(['./outputs/parameter_Morris/'], 'dir')
   mkdir(['./outputs/parameter_Morris/'])
end
if ~exist(['./outputs/parameter_sweep_ASM1seq/'], 'dir')
   mkdir(['./outputs/parameter_sweep_ASM1seq/'])
end
if ~exist(['./outputs/parameter_sweep_PBE_lambda/'], 'dir')
   mkdir(['./outputs/parameter_sweep_PBE_lambda/'])
end
if ~exist(['./outputs/plots/'], 'dir')
   mkdir(['./outputs/plots/'])
end
 
 switch P.IC_at
    case 'SS_computed'
        folderInput = ['outputs/data/Tsim_',num2str(P.T_sim),'/'];
        folderOutput = ['outputs/data/Tsim_',num2str(P.T_sim),'_SS/'];
    case 'data' 
        folderOutput = ['outputs/data/Tsim_',num2str(P.T_sim),'/'];
        folderInput = folderOutput;
    case 'startupScenario'
        folderInput = ['outputs/data/Tsim_',num2str(P.T_sim),'/'];
        folderOutput = ['outputs/data/Tsim_',num2str(P.T_sim),'_lowBiom/'];
 end
 
 % - Save data location
 if ~exist(folderInput, 'dir') && ~strcmp(IC_at,'data')
    error('First the steady state of the ASM1seq with the identified optimal parameter set needs to be computed. Rerun Main with IC_at = data and parSetup = chosenSetASM1seq')
 end
 if ~exist(folderOutput, 'dir')
   mkdir(folderOutput)
 end


 
 % - Plotting options 
 folderPlots = [folderOutput 'plots/'];
 if ~exist(folderPlots, 'dir')
   mkdir(folderPlots)
 end
 
 namesGlobal = who('global');  %# A cell array of variable names
 namesLocal = who();  
 namesGLocal = intersect(namesGlobal,namesLocal); 
 for iVar = 1:numel(namesGLocal)
  P.(namesGLocal{iVar}) = eval(namesGLocal{iVar});  % [EDITED]
 end
return

function [] = defineParams_WWT()
global N_tot V_0  spNames pres spLines  solubLines solubNames solubPlotted InitCondition
global nu_sosiaHET  mu_adh_ vMin_REAL vMax_REAL mu_miCol_v_REAL
global diamMin_T diamInitPart_T diamMax_T spColors InitVol_whole_aSp 
global HET AOB NOB sosiaHET EPS I mu_max NO2 NO3 S O2 NH4 K nu_HET Y S_in decay   
global HRT scaleDots   modesGrHET   speciesUnif nSpPBE nSpU
global speciesPBE nSp InitBiom_aSp  miColDensity sigma_adh_ stDev_v_wrt_mu_v 
global vMin_T vMax_T  switchTimes  dilInfl dilEffl dilPurge q_02_steady targetO2conc
global mu sigma mu_v_T sigma_v_T devSt_v_T  mu_adh_T  sigma_adh_T  devSt_adh_T  
global purgeWindow hoursOneDay minsOneHour Qinfl QEffl QPurge 

hoursOneDay = 24;
minsOneHour = 60;
%%% operational parameters
%% %%%%%%%%%%%%%%%%%%%%%%% WWT Parameters %%%%%%%%%%%%%%%%%%%%%% 

% indicate the minutes since the tControl = 0 at which we have a switch of phase
% phases . 1 . 2 3 4 . 2  3 . 4 2 . 3 . 4  5 .  6 
switchTimes = cumsum([50,15,25,117,15,25,117,15,25,118,140+3,55]); %data taken from excel
purgeWindow = 10;

% HRT = V_0/Q_cont
% SRT = 20 days
% Q1 = influx Q2=outflux_fluids Q3=outflux_solids // these are instantaneous fluxes
% Vd1 = volumeIn1d . V2=volOut1d_fluids . V3=volOut1d_solids
% Wd1 = window(fractions) of one day in which the respective flux is 
% solve system of equations Vd1 = Vd2+Vd3 -->Q1*W1=Q2*W2+Q3*W3
%                                      HRT*W2*Q2+ *W3*Q3 = 1000
% Q2 = (Q1*W1*SRT-V_0)/(W2*(SRT-HRT));
% Q3 = (V_0-Q1*W1*HRT)/(W3*(SRT-HRT));

V_0 = 1000; %@RealPar m^3 reactor volume % NB volume changes in phases 3 and 6
Q1 = 210;%@RealPar %m^3/d flux I have if I consider a continuous flux in time  
HRT = V_0/Q1; % 0.3; %days 
CycleTime = switchTimes(end);
tau = CycleTime/minsOneHour/hoursOneDay; %cycle time days/cycle
SRT = 20; % sludge retention time days %@RealPar

influx_tMins = 25*3; % Min/cycle
influx_tDays = influx_tMins/minsOneHour/hoursOneDay; %days/cycle
W1 = influx_tDays/tau; %[dimensionless] Fraction of time actually used for influx. 1/tau=2 is the number of cycles in one day

efflux_tMins = 55;% Min/cycle
efflux_tDays = efflux_tMins/minsOneHour/hoursOneDay; %days/cycle
W2 = efflux_tDays/tau; %Fraction of one day actually used for efflux. 1/tau=2 is the number of cycles in one day

purge_tMins = purgeWindow;% Min/cycle
purge_tDays = purge_tMins/minsOneHour/hoursOneDay; %days/cycle
W3 = purge_tDays/tau; %Fraction of one day actually used for efflux. 1/tau=2 is the number of cycles in one day

O2_tMins = 50+117+117+118;% Min/cycle
O2_tDays = O2_tMins/minsOneHour/hoursOneDay; %days/cycle
WO2 = O2_tDays/tau; %Fraction of one day actually used for efflux. 1/tau=2 is the number of cycles in one day
targetO2conc = 1.5/1000;% g/L == kg/m^3
q_02_steadymin = targetO2conc;% Kg/(m^3*(10mins))I want targetO2conc to be reached in dtO2~10 mins
dtO2 = 10; %mins
q_02_steady = q_02_steadymin*(1/(dtO2/minsOneHour/hoursOneDay)); % how many periods of 10 mins I have in one day
Q2 = V_0/HRT/W2;%1000m^3/HRT/(effective time in 1day)
Q3 = V_0/SRT/W3;%1000m^3/20days/(effective time in 1day)
dilInfl = Q1/V_0; % for the influx --> HRT/W1 * this is what I am doing combining the previous
dilEffl = Q2/V_0; % for the efflux --> HRT/W2
dilPurge = Q3/V_0; % for the purge --> SRT/W3
Qinfl = Q1;
QEffl = Q2;
QPurge = Q3; 
% solidsConc = 7; % g/L==kg/m^3
% solidDensity = 1060; % g/L==kg/m^3
% fracBiomassInSolids = 0.75;

% q_O2 = 0.5; %kg/(d*m^3)  **************missing************** 
% S_star_O2 = 0.009;%kg/m3  **************missing************** 


%growth parameters for indexes HET = 0,AOB=1,NOB=2. Taken from
mu_max = zeros(sosiaHET,1); %1/d converted below! 
mu_max(HET) = 6;
mu_max(sosiaHET) = 6;
mu_max(AOB) = 0.76;
mu_max(NOB) = 0.81;                             


%% %%%%%%%%%%%%%%%%%%%%%%% code settings %%%%%%%%%%%%%%%%%%%%%%
nSp = length(speciesModeled);
nSpPBE = length(speciesPBE);
nSpU = length(speciesUnif);
spNames = {'HET','AOB','NOB','sosiaHET'};
allBaseSpecies = [AOB HET NOB sosiaHET];
spLines = {'r-','c-','g-','r-'};
% spLines.log = {'r:','c:','g:','r:'};
C = linspecer(nSp); 
ll = 1;
for aSp = speciesModeled
 aSpN = char(spNames{aSp});
 spColors.(aSpN) = C(ll,:);
 ll = ll+1;
end
scaleDots = 'bo';
% scaleDots.log = 'y*';


%indexing for S vector
O2 = 1; S=2; NH4=3; NO2=4; NO3=5; EPS=6; I=7; % . EPS modeled here as a soluble as I am not interested to the PSD of EPS
solubles = [O2,S,NH4,NO2,NO3];
solubPlotted = solubles;
solubNames = {'O2','S','NH4','NO2','NO3','EPS','I'};
solubLines = {'-','-','-','-','-','-','-'};
% solubLines.log = {':',':',':',':',':',':',':'};
solPartic = [];
if modelEPS
 solPartic = EPS;
 if modelI
    solPartic = [solPartic, I]; % occhio se modello I do per scontato che sto modellizzando anche EPS
 end
end
solubPlotted = [solubles solPartic];
lastSub = solubPlotted(end);
solubNames = solubNames(solubPlotted);
pres = ismember(allBaseSpecies,speciesModeled);
modesGrHET = {'AER', 'NO2','NO3'};


decay = zeros(EPS,1); % 1/d converted below!
decay(HET) = 0.4; % 
decay(sosiaHET) = 0.4; 
decay(AOB) = 0.11;                                            
decay(NOB) = 0.11;                                            
% decay(EPS) = 0.17;

Y = zeros(I,1); % @@NT in principle Y should not need unit conversion
Y(HET) = 0.61; %gCOD/gCOD for HET 
Y(sosiaHET) = 0.61; %gCOD/gCOD for HET                                            
Y(AOB) = 0.33; %gCOD/gN for AOB and NOB                                            
Y(NOB) = 0.083; %gCOD/gN for AOB and NOB                                            
% Y(EPS) = 0.18; %gCOD/gN for AOB and NOB                                            

K = zeros(NO3,NOB); %kg/m^3 converted below!
K(O2,HET) = 0.81;
K(O2,sosiaHET) = 0.81;
K(O2,AOB) = 0.5e-3;
K(O2,NOB) = 0.68e-3;
K(S,HET) = 1e-2;
K(S,sosiaHET) = 1e-2;
K(NH4,AOB) = 1e-3;
K(NO2,HET) = 0.3e-3;
K(NO2,sosiaHET) = 0.3e-3;
K(NO2,NOB) = 1.3e-3;
K(NO3,HET) = 0.3e-3;
K(NO3,sosiaHET) = 0.3e-3;

nu_HET = 0.6; %dimensionless % ?Reduction factor in anoxic conditions
nu_sosiaHET = 0.6;%dimensionless % ?Reduction factor in anoxic conditions

S_in = zeros(lastSub,1); %?kg/m^3 %kg/m^3 converted below!
switch setup
 case 'AMBcurtis'   
    S_in(S) = 0.04;
    S_in(O2) = 0.005;
    S_in(NH4) = 0.04;
    S_in(NO2) = 0.001;
    S_in(NO3) = 0;
 case 'plantFlanders'
    S_in(S) = 1680; %kg/m^3
    S_in(O2) = 0; %  **************missing************** 
    S_in(NH4) = 0.105*1000; % kg/dm^3 --> kg/m^3
    S_in(NO2) = 0; %  **************missing************** 
    S_in(NO3) = 0;
end

if modelEPS
 Y(EPS) = 0.18; %gCOD/gN for AOB and NOB
 decay(EPS) = 0.17;
 S_in(EPS) = 0; % all'inizio metto ZERO o stessi grammi di EPS quanti grammi di HET, come nel modello?
end
if modelI
 Y(I) = 0.4;
 S_in(I) = 0;
end
if ~decayActive
 decay = zeros(EPS,1);
end

% setConversion();
% if convertUnits
% convertU();
% end

% DEFINE COORDINATES
%     _REAL_ : coordinates that appear in the plots
%     _T_ : coordinates calculated based on the option of what I want to represent with x2  [switch
%     InitCondition]
%     _N_ : normalized version 
vMin_REAL = 4/3*pi*(diamMin_T/2)^3; 
vMax_REAL = 4/3*pi*(diamMax_T/2)^3; 
mu_miCol_v_REAL = 4/3*pi*(diamInitPart_T/2)^3;
% mu_cellBiom_ = 1e-12; %g; https://hypertextbook.com/facts/2003/LouisSiu.shtml
% mu_cellVolume_ = (1/micronOneMeter)^3*0.4; % average cell volume in WW bacteria 0.4 microm^3 MEASURING BACTERIAL BIOMASS-COD IN WASTEWATER CONTAINING PARTICULATE MATTER

%% Define true _T values
% 1) for x1: mu_adh, dev_st_adh, sigma_adh
for aSp = speciesModeled 
 aSpN = char(spNames{aSp});
 mu_adh_T.(aSpN) = mu_adh_(aSp);
 sigma_adh_T.(aSpN) = sigma_adh_;
 devSt_adh_T.(aSpN) = sqrt(sigma_adh_T.(aSpN));
end

% 2) for x2: mu_v, dev_st_v, sigma_v (+ [vMin vMax], InitVol_whole_aSp, InitBiom_aSp, N_tot(initial number of particles))
for aSp = speciesPBE 
 aSpN = char(spNames{aSp});
 mu_miCol_Biom_aSp = mu_miCol_v_REAL*miColDensity;
 switch InitCondition
    case 'volume' 
        mu_v_T.(aSpN) = mu_miCol_v_REAL; % In this case they are all equal
        vMin_T.(aSpN) = vMin_REAL;
        vMax_T.(aSpN) = vMax_REAL;
        InitBiom_aSp = InitVol_whole_aSp*miColDensity;
        N_tot = ones(1,aSp)*(InitVol_whole_aSp/mu_miCol_v_REAL);
    case 'biomass' % impo!! in this case I am actually modeling microcolony biomass and not microcolony volume!
        mu_v_T.(aSpN) = mu_miCol_Biom_aSp;
        vMin_T.(aSpN) = vMin_REAL*miColDensity;
        vMax_T.(aSpN) = vMax_REAL*miColDensity;
        InitVol_whole_aSp = InitBiom_aSp/miColDensity;
        N_tot = ones(1,aSp)*(InitBiom_aSp/mu_miCol_Biom_aSp);
    case 'ntot' % impo!! here I impose N_tot and derive the InitBiom_aSp
        mu_v_T.(aSpN) = mu_miCol_Biom_aSp;
        vMin_T.(aSpN) = vMin_REAL*miColDensity;
        vMax_T.(aSpN) = vMax_REAL*miColDensity;
        InitBiom_aSp = N_tot(aSp)*mu_miCol_Biom_aSp;
        InitVol_whole_aSp = InitBiom_aSp/miColDensity;
    case 'ntot_save_InitBiom_aSp_sballa_miColDensity' % impo!! in this case I am actually modeling microcolony biomass and not microcolony volume!
        mu_miCol_Biom_aSp = InitBiom_aSp/N_tot(aSp);
        miColDensity = mu_miCol_Biom_aSp/mu_miCol_v_REAL;
        mu_v_T.(aSpN) = mu_miCol_Biom_aSp;
        vMin_T.(aSpN) = vMin_REAL*miColDensity;
        vMax_T.(aSpN) = vMax_REAL*miColDensity;
        InitVol_whole_aSp = InitBiom_aSp/miColDensity;
 end
 devSt_v_T.(aSpN) = mu_v_T.(aSpN)*stDev_v_wrt_mu_v;
 sigma_v_T.(aSpN) = devSt_v_T.(aSpN)^2;

% put both coordinate into singular mu and sigma matrixes
 mu.(aSpN) = [mu_adh_T.(aSpN) mu_v_T.(aSpN)]; % in cm  
 sigma.(aSpN) = [sigma_adh_T.(aSpN) sigma_v_T.(aSpN)].*eye(2); % in cm  !!NB .*eye(2) IS IMPORTANT, when I calculate the joint probability with mvnpdf I need also COVARIANCES
end
return

% function [] = setConversion()
% %%% Units conversions
% global mToMicrons daysTohrs kgToGr hoursOneDay grOneKg newOneMeter mmOneMeter micronOneMeter convertUnits
% mmOneMeter = 1e3;
% micronOneMeter = 1e6; 
% hoursOneDay = 24;
% grOneKg = 1000;
% LogicalStr = {'false', 'true'};
% if convertUnits
% if mToMicrons 
%     newOneMeter = micronOneMeter;
% else
% if mTomm
%     newOneMeter = mmOneMeter;
% end
% end
% disp(strcat(['NB units have been converted: m(meters) to microns [',LogicalStr{mToMicrons},'], to mm [',LogicalStr{mTomm},']; days to hrs [',LogicalStr{daysTohrs},']; kg to g [',LogicalStr{kgToGr},']']))
% end
% return

% function [] = convertU()
% global solPartic solubles q_O2 S_star_O2  diamMin_T diamInitPart_T diamMax_T 
% global mu_max K S_in decay      newOneMeter   HRT SDT   grOneKg miColDensity 
% %conversions
% miColDensity = miColDensity*(1/grOneKg)*(1e03)*(1e03); % from g/mL to Kg/m3:
% diamMin_T = diamMin_T  *newOneMeter;
% diamInitPart_T = diamInitPart_T  *newOneMeter;
% diamMax_T = diamMax_T  *newOneMeter; 
% % ------------
% q_O2 = q_O2*(grOneKg/(hoursOneDay*newOneMeter^3));   %g/(h*micron^3)      
% S_star_O2 = S_star_O2*(grOneKg/newOneMeter^3);  %g/micron^3      
% HRT = HRT*hoursOneDay;  %h                      
% SDT = HRT*(alpha+beta)/(beta*(1+alpha)); %h 
% S_in(solubles) = S_in(solubles) .*(grOneKg/newOneMeter^3); %g/micron^3 @@NT
% S_in(solPartic) = S_in(solPartic).*grOneKg;  %@@NT . NB different units in S_in!!
% K = K .*(grOneKg/newOneMeter^3); %g/micron^3  @@NT
% decay = decay .*(1/hoursOneDay); %1/h
% mu_max = mu_max .*(1/hoursOneDay); %1/h
% return


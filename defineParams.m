function [] = defineParams()
global convertUnits InitBiom_all kineticParams
global HET AOB NOB NO2 S O2 NH4 K S_in  S_0  

InitBiom_all = 5; % Kg
initProportions = [1,0.05,0.05]; % initial proportions among biomass components, wrt HET
initProportions = initProportions/sum(initProportions);
InitBiom_all = InitBiom_all*initProportions;
switch kineticParams
 case 'AMBcurtis'   
    [K_ABM,S_in_ABM] = setBaselineKineticParam();
    disp('remember that initial Biomass 0.02kg for each specie are set as initial biomass concentration, see paper')        
    disp(strcat(['Specific substrates:'])); 
    disp(strcat(['for HET is carbon conc = ',num2str(S_0(S)),', vs K(S,HET)=,',num2str(K(S,HET))]));      
    disp(strcat(['for AOB is NH4 conc = ',num2str(S_0(NH4)),', vs K(NH4,AOB)=,',num2str(K(NH4,AOB))]));        
    disp(strcat(['for NOB is NO2 conc = ',num2str(S_0(NO2)),', vs K(NO2,NOB)=,',num2str(K(NO2,NOB))]));        
    disp(strcat(['OXYGEN conc = ',num2str(S_0(O2)),': Kcat for HET=',num2str(K(O2,HET)),'; Kcat for AOB=',num2str(K(O2,AOB)),'; Kcat for NOB=',num2str(K(O2,NOB))]));
    S_in = S_in_ABM;
    S_0 = S_in;
    K = K_ABM;
 case {'actualSystem', 'ASM1'}
    S = initializeSmatrix();
    defineKineticParams();
    [X_in_ASM1,S_in_ASM1] = defineInletConcentrations();
    this_x_all_0_ASM1 = defineInitialConditions_t0(X_in_ASM1,S_in_ASM1);
end
setConversion();
if convertUnits
    convertU();
end
setMuSigma()
defineNtot()
return

function [] = setConversion()
% Units conversions
global mToMicrons gToKg daysTohrs kgToGr hoursOneDay grOneKg newOneMeter mmOneMeter micronOneMeter convertUnits minsOneHour
minsOneHour = 60;
mmOneMeter = 1e3;
micronOneMeter = 1e6; 
hoursOneDay = 24;
grOneKg = 1000;
gToKg = 1/grOneKg;
LogicalStr = {'false', 'true'};
if convertUnits
 if mToMicrons 
    newOneMeter = micronOneMeter;
 else
 if mTomm
    newOneMeter = mmOneMeter;
 end
 end
 disp(strcat(['NB units have been converted: m(meters) to microns [',LogicalStr{mToMicrons},'], to mm [',LogicalStr{mTomm},']; days to hrs [',LogicalStr{daysTohrs},']; kg to g [',LogicalStr{kgToGr},']']))
end
return

function [] = convertU()
    global solPartic solubles q_O2 S_star_O2 diamMin_T diamInitPart_T diamMax_T 
    global mu_max K S_in decay newOneMeter HRT SDT grOneKg miColDensity 
    
    % Conversions
    miColDensity = miColDensity * (1 / grOneKg) * (1e03) * (1e03); % from g/mL to Kg/m3
    diamMin_T = diamMin_T * newOneMeter;
    diamInitPart_T = diamInitPart_T * newOneMeter;
    diamMax_T = diamMax_T * newOneMeter; 
    
    % ------------
    q_O2 = q_O2 * (grOneKg / (hoursOneDay * newOneMeter^3)); % g/(h*micron^3)      
    S_star_O2 = S_star_O2 * (grOneKg / newOneMeter^3); % g/micron^3      
    HRT = HRT * hoursOneDay; % h                      
    SDT = HRT * (alpha + beta) / (beta * (1 + alpha)); % h 
    S_in(solubles) = S_in(solubles) .* (grOneKg / newOneMeter^3); % g/micron^3 
    S_in(solPartic) = S_in(solPartic) .* grOneKg; % NB different units in S_in!
    K = K .* (grOneKg / newOneMeter^3); % g/micron^3 
    decay = decay .* (1 / hoursOneDay); % 1/h
    mu_max = mu_max .* (1 / hoursOneDay); % 1/h    
return


function [K_ABM,S_in_ABM] = setBaselineKineticParam()
global Y decay decayActive nu_HET mu_max HET AOB NOB O2 S NH4 NO2 NO3 I EPS lastSub modelEPS modelI
%growth parameters for indexes HET = 0,AOB=1,NOB=2. Taken from
mu_max = zeros(NOB,1); %1/d converted below! 
mu_max(HET) = 6;
mu_max(AOB) = 0.76;
mu_max(NOB) = 0.81; 

decay = zeros(EPS,1); % 1/d converted below!
decay(HET) = 0.4; % 
decay(AOB) = 0.11;                                            
decay(NOB) = 0.11;                                            

Y = zeros(I,1); % Y does not need unit conversion
Y(HET) = 0.61; % gCOD/gCOD for HET 
Y(AOB) = 0.33; % gCOD/gN for AOB and NOB                                            
Y(NOB) = 0.083; % gCOD/gN for AOB and NOB                                            
nu_HET = 0.6; % dimensionless Reduction factor in anoxic conditions

K_ABM = zeros(NO3,NOB); % kg/m^3 
K_ABM(O2,HET) = 0.81; 
K_ABM(O2,AOB) = 0.5e-3;
K_ABM(O2,NOB) = 0.68e-3;
K_ABM(S,HET) = 1e-2;
K_ABM(NH4,AOB) = 1e-3;
K_ABM(NO2,HET) = 0.3e-3;
K_ABM(NO2,NOB) = 1.3e-3;
K_ABM(NO3,HET) = 0.3e-3;

% SUBSTRATES
S_in_ABM = zeros(lastSub,1); %?kg/m^3 %kg/m^3 converted below!
S_in_ABM(O2) = 0.005;
S_in_ABM(S) = 0.04;
S_in_ABM(NH4) = 0.04;
S_in_ABM(NO2) = 0.001;
S_in_ABM(NO3) = 0;

if modelEPS
 Y(EPS) = 0.18; % gCOD/gN for AOB and NOB
 decay(EPS) = 0.17;
 S_in_ABM(EPS) = 0; 
end
if modelI
 Y(I) = 0.4;
 S_in_ABM(I) = 0;
end
if ~decayActive
 decay = zeros(EPS,1);
end
return

function [] = setMuSigma()
global mu mu_adh_T sigma_adh_T devSt_adh_T mu_v_T diamInitPart_T diamMax_T InitVol_all InitBiom_all
global vMin_T vMax_T ModelBiomass cells1KgBiom speciesPBE stDev_v_wrt_mu_v sigma_v_T sigma_adh_ biomassDensity
global sigma devSt_v_T mu_adh_ spNames   diamMin_T 
% DEFINE COORDINATES
%     _T_ : coordinates expressed as volume or biomass
%     _N_ : normalized version 

%% Define true _T values
for aSp = speciesPBE 
 aSpN = char(spNames{aSp});
 
% 1) for x1: mu_adh, dev_st_adh, sigma_adh
 mu_adh_T.(aSpN) = mu_adh_(aSp);
 sigma_adh_T.(aSpN) = sigma_adh_;
 devSt_adh_T.(aSpN) = sqrt(sigma_adh_T.(aSpN));

% 2) for x2: mu_v, dev_st_v, sigma_v (+ [vMin vMax], InitVol_whole_aSp, InitBiom_aSp, N_tot(initial number of particles))% Particle volume coordinate
 diamMin_T = 0; % lower and upper particle volumes
 diamInitPart_T = 1e-6; % m
 diamMax_T = 200e-6; %150e-06; % m
 stDev_v_wrt_mu_v = 1/100;% gives me the proportion of stdev/mu
 cells1KgBiom = 5E14;%cells per kg of biomass, from cells1g=5E11;% 
% biomassDensity = 1.060; % g/L==kg/m^3

 mu_v_T.(aSpN) = 4/3*pi*(diamInitPart_T/2)^3; % In this case they are all equal
 vMin_T.(aSpN) = 4/3*pi*(diamMin_T/2)^3;
 vMax_T.(aSpN) = 4/3*pi*(diamMax_T/2)^3; 
 
 biomassDensity = (1/cells1KgBiom)/mu_v_T.(aSpN); %Kg/m^3 density of a single cell
 InitVol_all = InitBiom_all/biomassDensity;
 disp(strcat(['Please note that with these parameters the cell density results = ',num2str(biomassDensity)]));
 
% singleCellBiomass = mu_v_T.(aSpN)*biomassDensity;
 
 if ModelBiomass % impo!! in this case I am actually modeling microcolony biomass and not microcolony volume!
    mu_v_T.(aSpN) = mu_v_T.(aSpN)*miColDensity;
    vMin_T.(aSpN) = vMin_T.(aSpN)*miColDensity;
    vMax_T.(aSpN) = vMax_T.(aSpN)*miColDensity;
 end
 devSt_v_T.(aSpN) = mu_v_T.(aSpN)*stDev_v_wrt_mu_v;
 sigma_v_T.(aSpN) = devSt_v_T.(aSpN)^2;

% put both coordinate into singular mu and sigma matrixes
 mu.(aSpN) = [mu_adh_T.(aSpN) mu_v_T.(aSpN)]; % in cm  
 sigma.(aSpN) = [sigma_adh_T.(aSpN) sigma_v_T.(aSpN)].*eye(2); % in cm  !!NB .*eye(2) IS IMPORTANT, when I calculate the joint probability with mvnpdf I need also COVARIANCES
end
return

function [] = defineNtot()
global N_tot_all InitBiom_all cells1KgBiom 
 N_tot_all = InitBiom_all*cells1KgBiom;
return




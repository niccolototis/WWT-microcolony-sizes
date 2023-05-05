function [P] = ASM1seqParameters(P)
% P = defineKineticParams(P);
[X_c_in,S_c_in] = defineInletConcentrations(P);
P.X_c_in = X_c_in;
P.S_c_in = S_c_in;
% P = defineInitialConditions_t0(P);
% %-- In case I have inflow of particles that are modeled with the PBE
% %I have to take that into account later in the ode_fun of the PBE, to add
% %particle to the system during the feed phase
% partIn = X_c_in(P.speciesPBE);
% P.particleFlowIn = false;
% if length(partIn(partIn>0))>0
% P.particleFlowIn = true;
% end

% namesGlobal = who('global'); %# A cell array of variable names
% namesLocal = who(); 
% namesGLocal = intersect(namesGlobal,namesLocal); 
% for iVar = 1:numel(namesGLocal)
% P.(namesGLocal{iVar}) = eval(namesGLocal{iVar}); % [EDITED]
% end
return

function P = defineKineticParams(P)
 P.MuH = 6.0; % per day; heterotroph growth -  alternative 3.0
 P.Ks = 20; %0.6; %20.0; @ASM1 % g COD/m^3; COD half-saturation coefficient - 0.6 - 2
 P.KOH = 0.20; % g O2/m^3; oxygen switch  
 P.KNO = 0.50; % g NO3-N/m^3; nitrate switch / 4
 P.bH = 0.62;%0.2; %0.62 @ASM1; % per day; heterotroph decay / 0.4 / 0.2 (aerobic) or 0.1 (anoxic) / 0.3 (aerobic) or 0.1 (anoxic, batch) or 0.15 (anoxic, full-scale)	/ 0.3	
 P.bA = 0.10; % per day; autotroph decay  / 0.28 (aerobic test); 0.07 (anoxic test)	/ 0.15
 P.Etag = 0.8; % anoxic growth discount /0.4 /0.6
 P.Etah = 0.4; % anoxic hydrolysis discount 
 P.kh = 3.0; % g slowly biodegradable COD/g cell COD/day; hydrolysis rate
 P.KX = 0.03; % g slowly biodegradable COD/g cell COD; hydrolysis half-saturation coeff.
 P.MuA = 0.8;%2; %0.80 @ASM1; % per day; autotroph growth - 0.76 / 1.05 - 1.4  / 0.8-1.1  / 0.4	-  2.16
 P.KNH = 1.0; % g NH3-N/m^3; ammonia half-saturation coefficient  0.01 - 0.5
 P.KOA = 0.4; % g O2/m^3; oxygen switch / 0.5
 P.ka = 0.08; % m^3/g COD/d; ammonification rate 
return

function [X_c_in,S_c_in] = defineInletConcentrations(P) 
X_c_in = zeros(P.nPart,1); %/ mg COD L-1
S_c_in = zeros(P.nSolub,1); %/ mg COD L-1
switch P.kineticParams 
 case 'ASM1_original'
    %    g/m3 = mg/L
    X_c_in(P.HET) = 30;  %Active heterotrophic biomass, Xb,h /  mg COD L-1
    X_c_in(P.AUT) = 1;  %Active autotrophic biomass, Xb,a /  mg COD L-1
    X_c_in(P.xI) = 25;  %Particulate inert organic matter, Xi /  mg COD L-1
    X_c_in(P.xS) = 125;  %Slowly biodegradable substrate, Xs /  mg COD L-1
    X_c_in(P.xP) = 0;  % Particulate decay products, Xp /  mg COD L-1
    X_c_in(P.xND) = 0;  %Particulate biodegradable organic nitrogen, Xnd /  mg COD L-1
    %    g/m3 = mg/L
    S_c_in(P.sI) = 30; % Soluble inert organic matter, Si /  mg COD L-1
    S_c_in(P.sS) = 30; %Readily biodegradable substrate, Ss /  mg COD L-1
    S_c_in(P.sNH4) = 16; %Ammonia-nitrogen, Snh /  mg COD L-1
    S_c_in(P.sNO) = 0; %Nitrate- and nitrite- nitrogen, Sno /  mg COD L-1
    S_c_in(P.sND) = 0; %Soluble biodegradable organic nitrogen, Snd /  mg COD L-1
    S_c_in(P.sALK) = 5; % Alkalinity - molar units, Salk /  mg COD L-1
    S_c_in(P.sO2) = 0; 
 case 'actualSystem'
    %    g/m3 = mg/L
    X_c_in(P.HET) = 0;  %Active heterotrophic biomass, Xb,h /  mg COD L-1
    X_c_in(P.AUT) = 0;  %Active autotrophic biomass, Xb,a /  mg COD L-1
    X_c_in(P.xI) = 25;  %Particulate inert organic matter, Xi /  mg COD L-1
    X_c_in(P.xS) = 125;  %Slowly biodegradable substrate, Xs /  mg COD L-1
    X_c_in(P.xP) = 0;  % Particulate decay products, Xp /  mg COD L-1
    X_c_in(P.xND) = 0;  %Particulate biodegradable organic nitrogen, Xnd /  mg COD L-1   
    %    g/m3 = mg/L
    S_c_in(P.sI) = 30; % Soluble inert organic matter, Si /  mg COD L-1
    S_c_in(P.sS) = 1680; %Readily biodegradable substrate, Ss /  mg COD L-1
    S_c_in(P.sNH4) = 105; %Ammonia-nitrogen, Snh /  mg COD L-1
    S_c_in(P.sNO) = 0; %Nitrate- and nitrite- nitrogen, Sno /  mg COD L-1
    S_c_in(P.sND) = 0; %Soluble biodegradable organic nitrogen, Snd /  mg COD L-1
%     S_c_in(P.sALK) = 5; % Alkalinity - molar units, Salk /  mg COD L-1
    S_c_in(P.sO2) = 0; 
end
if P.variableInGrams
 % do nothing, these remain concentrations, and are used in the odefun
 % as such
end
return

function P = defineInitialConditions_t0(P)
global HET AUT xI xS xP xND sI sS sO2 sNH4 sNO sND sALK targetO2conc nPart nSolub V_0 kineticParams
global variableInGrams InitBiom_all 
 X_0_ASM1 = zeros(P.nPart,1);
 S_0_ASM1 = zeros(P.nSolub,1);
 
 switch P.kineticParams 
 case {'ASM1_original', 'actualSystem'}
    %    g/m3 = mg/L
    X_0_ASM1(P.HET) = 3969;  %Active heterotrophic biomass, Xb,h
    X_0_ASM1(P.AUT) = 441;  %Active autotrophic biomass, Xb,a
    X_0_ASM1(P.xI) = P.X_c_in(P.xI);  %Particulate inert organic matter, Xi
    X_0_ASM1(P.xS) = P.X_c_in(P.xS);  %Slowly biodegradable substrate, Xs
    X_0_ASM1(P.xP) = P.X_c_in(P.xP);  % Particulate decay products, Xp
    X_0_ASM1(P.xND) = P.X_c_in(P.xND);  %Particulate biodegradable organic nitrogen, Xnd
    %    g/m3 = mg/L
    S_0_ASM1(P.sI) = P.S_c_in(P.sI); % Soluble inert organic matter, Si
    S_0_ASM1(P.sS) = 10;  %Readily biodegradable substrate, Ss
    S_0_ASM1(P.sNH4) = 16;  %Ammonia-nitrogen, Snh
    S_0_ASM1(P.sNO) = 1;  %Nitrate- and nitrite- nitrogen, Sno
    S_0_ASM1(P.sND) = 1;  %Soluble biodegradable organic nitrogen, Snd
%     S_0_ASM1(P.sALK) = 5; % Alkalinity - molar units, Salk
    S_0_ASM1(P.sO2) = P.O2_c_ref;%targetO2conc; 
 end 
 %from concentrations to weight
 X_0_ASM1_weight = X_0_ASM1*P.V_0; % grams
 S_0_ASM1_weight = S_0_ASM1*P.V_0;
 P.InitBiom_all = X_0_ASM1_weight([P.HET P.AUT]); % this is needed for the PBE, needs to be in weight, not concentrations
 if P.variableInGrams % from concentrations to grams
    X_0_ASM1 = X_0_ASM1_weight;
    S_0_ASM1 = S_0_ASM1_weight;
 end
 this_x_all_0 = [X_0_ASM1; S_0_ASM1;P.V_0];
 P.this_x_all_0 = this_x_all_0;
return


% function [] = temperatureCorrections()
% global mu_max decay K Y nu_HET HET NOB 
% 
% T_plant = 24;
% theta = [1.07,1.1,1.1];
% tCorr = exp(theta.*(T_plant-20)); % correction factor for temp
% 
% mu_max = mu_max.*tCorr; % apply temp correction
% decay(HET:NOB) = decay(HET:NOB).*tCorr; % apply temp correction
% K = K.*tCorr; % apply temp correction
% nu_HET = nu_HET*tCorr(HET);
% return






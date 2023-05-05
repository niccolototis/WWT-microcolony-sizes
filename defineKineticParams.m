function P = defineKineticParams(P)
% [1] Ghorbani, M., and C. Eskicioglu. "Application of the International Water Association activated sludge models to describe aerobic sludge digestion." Environmental technology 32.16 (2011): 1923-1938.
% select the scaled parameter change delta \in [0,1]

% - these are to be used in defineInletConcentrations
% Source: ﻿Updated Activated Sludge Model n°1 Parameter Values for Improved Prediction of Nitrogen Removal in Activated Sludge Processes: Validation at 13 Full-scale Plants
% Author(s): Jean-Marc Choubert, Anne-Emmanuelle Stricker, Aurélien Marquot, Yvan Racault, Sylvie Gillot and Alain Héduit
% COD gO2/m3 = 500(+/-180) ; 634(+/-315)
% COD fractionation
% sS = 20%
% xS = 59%
% sI = 4%
% xI = 17%
P.f_Xi_in = 17;
P.f_Xs_in = 59;
P.f_Si_in = 4;
P.f_Ss_in = 20;

% - These parameters are used to set the initial conditions
P.biom_conc = 7000;% [g MLSS/m^3] biomass concentration in reactor at the SS @actualSystem
P.active_biom_in_gCOD = 0.3; % [g/g] active biomass_active / biomass_total biomass @assumed
P.gHET_in_total = 0.98;% [g/g]
P.gCOD_in_gVSS = 1.4;% [g COD/g VSS] @assumed based on Contreras 2012	
P.gVSS_in_gMLSS = 0.75;% [gVSS/g MLSS]

switch P.modelKinPars 
 case {'actualSystem'} % these are the original parameters of the ASM1 model at 20 deg Celsius        
    P.YA = 0.24; % g cell COD formed/ g N oxidized
    P.YH = 0.67; % g cell COD formed/ g COD oxidized
    
    P.fp = 0.08; % dimensionless, fraction of biomass going to inert products 
    P.iXB = 0.086; % g N in biomass / g COD in biomass
    P.iXP = 0.06; % g N in endogenous mass/ g COD in endogenous mass

    P.MuH = 6.0; % per day; heterotroph growth
    P.bH = 0.62; % per day; heterotroph decay
    P.Ks = 20.0; % g COD/m^3; COD half-saturation coefficient
    P.KOH = 0.20; % g O2/m^3; oxygen switch
    P.KNO = 0.50; % g NO3-N/m^3; nitrate switch    
    
    P.Etag = 0.8; % anoxic growth discount
    P.Etah = 0.4; % anoxic hydrolysis discount
    
    % hydrolysis of trapped particulate substrate
    P.kh = 3.0; % g slowly biodegradable COD/g cell COD/day; hydrolysis rate
    P.KX = 0.03; % g slowly biodegradable COD/g cell COD; hydrolysis half-saturation coeff.
    
    P.MuA = 0.80; % per day; autotroph growth
    P.bA = 0.10; % per day; autotroph decay
    % KOH = 0.20; % g O2/m^3; oxygen switch
    P.KNH = 1.0; % g NH3-N/m^3; ammonia half-saturation coefficient
    
    % hydrolysis of trapped particulate nitrogen
    P.KOA = 0.4; % g O2/m^3; oxygen switch
    P.ka = 0.08; % m^3/g COD/d; ammonification rate

    %@NB KNH_mod is introduced to make the growth rate of HET dependent on NH, otherwise it would be depleted below zero
    P.KNH_mod = 1e-6;       

 case {'ASM1_original'}
 % Default values from Table 5, p. 24
    P.YA = 0.24; % g cell COD formed/ g N oxidized
    P.YH = 0.67; % g cell COD formed/ g COD oxidized
    P.fp = 0.08; % dimensionless, fraction of biomass going to inert products 
    P.iXB = 0.086; % g N in biomass / g COD in biomass
    P.iXP = 0.06; % g N in endogenous mass/ g COD in endogenous mass

    P.MuH = 6.0; % per day; heterotroph growth
    P.bH = 0.62; % per day; heterotroph decay
    P.Ks = 20.0; % g COD/m^3; COD half-saturation coefficient
    P.KOH = 0.20; % g O2/m^3; oxygen switch
    P.KNO = 0.50; % g NO3-N/m^3; nitrate switch    
    
    P.Etag = 0.8; % anoxic growth discount
    P.Etah = 0.4; % anoxic hydrolysis discount
    
    % hydrolysis of trapped particulate substrate
    P.kh = 3.0; % g slowly biodegradable COD/g cell COD/day; hydrolysis rate
    P.KX = 0.03; % g slowly biodegradable COD/g cell COD; hydrolysis half-saturation coeff.
    
    P.MuA = 0.80; % per day; autotroph growth
    P.bA = 0.10; % per day; autotroph decay
    % KOH = 0.20; % g O2/m^3; oxygen switch
    P.KNH = 1.0; % g NH3-N/m^3; ammonia half-saturation coefficient
    
    % hydrolysis of trapped particulate nitrogen
    P.KOA = 0.4; % g O2/m^3; oxygen switch
    P.ka = 0.08; % m^3/g COD/d; ammonification rate

end
return

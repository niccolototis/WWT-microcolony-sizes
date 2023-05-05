function P = defineInletConcentrations(P)
% To be noted: if variableInGrams == true do nothing, these remain
% concentrations, and are used in the odefun as such

% Source: ﻿Updated Activated Sludge Model n°1
% Parameter Values for Improved Prediction of Nitrogen Removal in Activated
% Sludge Processes: Validation at 13 Full-scale Plants Author(s): Jean-Marc
% Choubert, Anne-Emmanuelle Stricker, Aurélien Marquot, Yvan Racault,
% Sylvie Gillot and Alain Héduit

% COD gO2/m3 = 500(+/-180) ; 634(+/-315)
% COD fractionation
% sS = 20%
% xS = 59%
% sI = 4%
% xI = 17%
% BOD gO2/m3 = 220(+/-100) ; 302(+/-170)
% TKN gN/m3 = 43(+/-14) ; 52(+/-23)
% NH4-N gN/m3 = 28(+/-12) ; nd
% TSS gTSS/m3 = 236(+/-123) ; 302(+/-170)
% VSS % = 78(+/-5) ; nd

P.X_c_in = zeros(P.nPart,1); % [mg COD/L]
P.S_c_in = zeros(P.nSolub,1); % [mg COD/L]
 
switch P.modelKinPars 
 case 'actualSystem'
    %    [g/m3 = mg/L]
    P.X_c_in(P.HET) = 0;  %Active heterotrophic biomass, Xb,h /  [mg COD/L]
    P.X_c_in(P.AUT) = 0;  %Active autotrophic biomass, Xb,a /  [mg COD/L]
    P.X_c_in(P.xI) = 1680*(P.f_Xi_in/100);%25;  %Particulate inert organic matter, Xi /  [mg COD/L]
    P.X_c_in(P.xS) = 1680*(P.f_Xs_in/100);%125;  %Slowly biodegradable substrate, Xs /  [mg COD/L]   
    P.X_c_in(P.xP) = 0;  % Particulate decay products, Xp /  [mg COD/L]
    P.X_c_in(P.xND) = 0;  %Particulate biodegradable organic nitrogen, Xnd /  [mg COD/L]   
    %    g/m3 = mg/L
    P.S_c_in(P.sI) = 1680*(P.f_Si_in/100);%30; % Soluble inert organic matter, Si /  [mg COD/L]
    P.S_c_in(P.sS) = 1680*(P.f_Ss_in/100); %Readily biodegradable substrate, Ss /  [mg COD/L]
    P.S_c_in(P.sNH4) = 105; %Ammonia-nitrogen, Snh /  [mg COD/L]
    P.S_c_in(P.sNO) = 0; %Nitrate- and nitrite- nitrogen, Sno /  [mg COD/L]
    P.S_c_in(P.sND) = 0; %Soluble biodegradable organic nitrogen, Snd /  [mg COD/L]
    P.S_c_in(P.sO2) = 0; 
    
 case 'ASM1_original'
    %    [g/m3 = mg/L]
    P.X_c_in(P.HET) = 30;  %Active heterotrophic biomass, Xb,h /  [mg COD/L]
    P.X_c_in(P.AUT) = 1;  %Active autotrophic biomass, Xb,a /  [mg COD/L]
    P.X_c_in(P.xI) = 25;  %Particulate inert organic matter, Xi /  [mg COD/L]
    P.X_c_in(P.xS) = 125;  %Slowly biodegradable substrate, Xs /  [mg COD/L]
    P.X_c_in(P.xP) = 0;  % Particulate decay products, Xp /  [mg COD/L]
    P.X_c_in(P.xND) = 0;  %Particulate biodegradable organic nitrogen, Xnd /  [mg COD/L]
    %    g/m3 = mg/L
    P.S_c_in(P.sI) = 30; % Soluble inert organic matter, Si /  [mg COD/L]
    P.S_c_in(P.sS) = 30; %Readily biodegradable substrate, Ss /  [mg COD/L]
    P.S_c_in(P.sNH4) = 16; %Ammonia-nitrogen, Snh /  [mg COD/L]
    P.S_c_in(P.sNO) = 0; %Nitrate- and nitrite- nitrogen, Sno /  [mg COD/L]
    P.S_c_in(P.sND) = 0; %Soluble biodegradable organic nitrogen, Snd /  [mg COD/L]
    P.S_c_in(P.sO2) = 0; 
end
return




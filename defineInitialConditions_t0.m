function P = defineInitialConditions_t0(P,varargin)
% Note that initial conditions are defined either as concentrations or
% grams depending on P.variableInGrams
global InitBiom_all folderInput
switch P.parSetup
 case 'chosenSetASM1seq'
    label = '';
 case 'baseline'
    label = 'BASE_';
 case 'EEMorris'
    label = 'BASE_';
 case 'paramSweep'
    label = 'BASE_';
 otherwise
    ME = MException('either one of chosenSetASM1seq or baseline need to be selected');
    throw(ME);
end 
switch P.IC_at
 case 'SS_computed'
    switch P.parSetup
        case 'chosenSetASM1seq'
            % - Load previously saved data
            try
                load([folderInput,'resASM1seq_',label,num2str(P.T_sim),'.mat'],'X_c_res','S_c_res');
            catch ME
                error('If the input data are not available probably need to run the model with IC_at = data option')
            end
            X_c_0 = X_c_res(end,:)';
            S_c_0 = S_c_res(end,:)';
    end
    
 case 'data' %{'data','SS_computed','startupScenario'};
    X_c_0 = zeros(P.nPart,1);
    S_c_0 = zeros(P.nSolub,1);
    X_c_0(P.HET) = P.biom_conc*P.active_biom_in_gCOD*P.gHET_in_total*P.gCOD_in_gVSS*P.gVSS_in_gMLSS; % Active heterotrophic biomass, Xb,h
    X_c_0(P.AUT) = P.biom_conc*P.active_biom_in_gCOD*(1-P.gHET_in_total)*P.gCOD_in_gVSS*P.gVSS_in_gMLSS; %22; %Active autotrophic biomass, Xb,a
    X_c_0(P.xI) = P.X_c_in(P.xI); %Particulate inert organic matter, Xi
    X_c_0(P.xS) = P.X_c_in(P.xS); %Slowly biodegradable substrate, Xs
    X_c_0(P.xP) = P.X_c_in(P.xP); % Particulate decay products, Xp
    X_c_0(P.xND) = P.X_c_in(P.xND); %Particulate biodegradable organic nitrogen, Xnd
    %    [g/m^3 = mg/L]
    S_c_0(P.sI) = P.S_c_in(P.sI); % Soluble inert organic matter, Si
    S_c_0(P.sS) = 10; %Readily biodegradable substrate, Ss
    S_c_0(P.sNH4) = 16; %Ammonia-nitrogen, Snh
    S_c_0(P.sNO) = 1; %Nitrate- and nitrite- nitrogen, Sno
    S_c_0(P.sND) = 1; %Soluble biodegradable organic nitrogen, Snd
    %     S_c_0(P.sALK) = 5; % Alkalinity - molar units, Salk
    S_c_0(P.sO2) = P.O2_c_ref; % defined in setOperationParams.m
 case 'startupScenario'
    load([folderInput,'resASM1seq_',label,num2str(P.T_sim),'.mat'],'X_c_res','S_c_res');
    X_c_0 = X_c_res(end,:)';
    S_c_0 = S_c_res(end,:)';
    X_c_0(P.HET) = P.biomFract0*X_c_0(P.HET);
    X_c_0(P.AUT) = P.biomFract0*X_c_0(P.AUT);
end

%% from concentrations to grams
X_g_0 = X_c_0*P.V_0; % [g]
S_g_0 = S_c_0*P.V_0;
InitBiom_all = X_g_0([P.HET P.AUT]); % this is needed for the PBE, needs to be in g, not concentrations
P.InitBiom_all = InitBiom_all;
if P.variableInGrams 
 X_0 = X_g_0;
 S_0 = S_g_0;
else
 X_0 = X_c_0;
 S_0 = S_c_0;
end
P.this_x_all_0 = [X_0;S_0;P.V_0];
return


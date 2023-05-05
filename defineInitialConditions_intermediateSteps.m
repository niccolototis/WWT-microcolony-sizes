function P=defineInitialConditions_intermediateSteps(P)
X=P.this_x_all_0(1:P.nPart);
S=P.this_x_all_0((P.nPart+1):(P.nPart+P.nSolub));
V=P.this_x_all_0(end);  

% - during the aerobic phases Oxigen is kept at a reference concentration
if P.onAer 
    if P.variableInGrams
        S(P.sO2)=P.O2_c_ref*V;
    else
        S(P.sO2)=P.O2_c_ref;
    end
end

% - compute the flow rates and volume changes during onPurge and onEfflux
if P.washOutParticles
if P.onPurge
if P.variableInGrams
    X_g=X;
else
    % X has been defined as concentrations
    X_g=X*V;                   
end    

% - The amount of sludge purged per cycle is computed knowing that in SRT
% days the whole sludge would be purged
% The volume depletion is obtained from the sludge grams and
% its density    
X_g_pu_NextCy=X_g*(1/P.SRT)*(1/P.cyclesOneDay); % [g/cycle]
X_g_tot_pu_NextCy=sum(X_g_pu_NextCy);
P.V_pu_NextCy=X_g_tot_pu_NextCy/P.sludgeDensity;
P.V_ef_NextCy=P.V_in_1Cy-P.V_pu_NextCy;
QEffl=P.V_ef_NextCy/P.efflux_tDays;            
P.dVdt_ef=-QEffl;
if P.instantSludgeRemoval  
    % particulate components are purged instantanously
    X_g=X_g-X_g_pu_NextCy; 
    V=V-P.V_pu_NextCy; 
    if P.variableInGrams
        X=X_g;
    else
        X_c=X_g/V;
        X=X_c;
    end               
else
    % particulate components are purged continuously during the purge phase
    % dX_g_pu QPurge are computed before the purge phase starts
    P.dX_g_pu=-(X_g/P.SRT)/P.purge_tDays; 
    QPurge=P.V_pu_NextCy/P.purge_tDays;
    P.dVdt_pu=-QPurge;
end
end
else % In case I am considering I have no purge
    P.V_pu_NextCy=0;
    if P.influxActive
        P.V_ef_NextCy=P.V_in_1Cy-P.V_pu_NextCy;
        QEffl=P.V_ef_NextCy/P.efflux_tDays;
        P.dVdt_ef=-QEffl; %QEffl=Q2
    end
end
P.this_x_all_0=[X;S;V];
return

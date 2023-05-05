function [dxdt] = odeFunASM1seq(t, x_all,P) 
% The ODE system is written coherently to how it is presented in the
% original ASM1 publication [Henze, 2000]
global step steP
steP = steP+1; 
step = step+1;

x_all = fixLowerThreshold_fast(x_all,P); 

% %particulate components: HET AUT xI xS xP xND
X = x_all(1:P.nPart);

% %soluble components: sI sS sO2 sNH4 sNO sND % sALK leaving ALK out
S = x_all((P.nPart+1):(P.nPart+P.nSolub));
V = x_all(end);
if P.variableInGrams
 S_c = S/V; % converting S_ctrates to concentrations to compute the Monod kinetics
 S_g = S;
 X_g = X;
else
 S_c = S;
 S_g = S*V;
 X_g = X*V;
end

dxdt = x_all*0;
dXdt = X*0;
dSdt = S_c*0;
dVdt = 0;

%%
if P.implementGrowth
 Rho = zeros(8,1);
 
 % Calculate Rho == growth rate
	Rho(1) = P.MuH*(S_c(P.sS)/(P.Ks+S_c(P.sS)))*(S_c(P.sO2)/(P.KOH+S_c(P.sO2)))*X(P.HET);
 Rho(1) = Rho(1)*(S_c(P.sNH4)/(P.KNH_mod+S_c(P.sNH4))); % NB modified wrt ASM1 model! Introduced the dependency from sNH4, which otherwise would be depleted beyond zero 
	
 Rho(2) = P.MuH*(S_c(P.sS)/(P.Ks+S_c(P.sS)))*(P.KOH/(P.KOH+S_c(P.sO2)))*(S_c(P.sNO)/(P.KNO+S_c(P.sNO)))*P.Etag*X(P.HET);    
 Rho(2) = Rho(2)*(S_c(P.sNH4)/(P.KNH_mod+S_c(P.sNH4))); % NB modified wrt ASM1 model! Introduced the dependency from sNH4, which otherwise would be depleted beyond zero 

 Rho(3) = P.MuA*(S_c(P.sNH4)/(P.KNH+S_c(P.sNH4)))*(S_c(P.sO2)/(P.KOH+S_c(P.sO2)))*X(P.AUT);
	
 Rho(4) = P.bH*X(P.HET); 
	Rho(5) = P.bA*X(P.AUT);
	Rho(6) = P.ka*S_c(P.sND)*X(P.HET);
	Rho(7) = P.kh*((X(P.xS)/X(P.HET))/(P.KX+(X(P.xS)/X(P.HET))))*((S_c(P.sO2)/(P.KOH+S_c(P.sO2)))+P.Etah*(P.KOH/(P.KOH+S_c(P.sO2)))*(S_c(P.sNO)/(P.KNO+S_c(P.sNO))))*X(P.HET);
	Rho(8) = Rho(7)*X(P.xND)/X(P.xS);
 
	r = P.Smat'*Rho; 
 
 
 dXdt(P.HET) = r(5);
 dXdt(P.AUT) = r(6);
 dXdt(P.xI) = r(3);
 dXdt(P.xS) = r(4);
 dXdt(P.xP) = r(7);
 dXdt(P.xND) = r(12);
 
 dSdt(P.sI) = r(1);
 dSdt(P.sS) = r(2);
 dSdt(P.sO2) = r(8);
 dSdt(P.sNH4) = r(10);
 dSdt(P.sNO) = r(9);
 dSdt(P.sND) = r(11);
 % dSdt(sALK) = r(13); leaving alk out
end

%% Operational cycles control 
if P.onFeed && P.influxActive
 dVdt = P.QInfl;
 if P.variableInGrams % NB both P.S_c_in and P.X_c_in are concentrations
    dSdt = dSdt + P.S_c_in*P.QInfl; 
    dXdt = dXdt + P.X_c_in*P.QInfl;
 else
    P.dilInfl = P.QInfl/V;
    dSdt = dSdt + P.S_c_in*P.dilInfl - S_c*P.dilInfl; 
    dXdt = dXdt + P.X_c_in*P.dilInfl - X*P.dilInfl;
 end
end

if P.washOutParticles
 if P.onPurge % washing out just I and EPS    
    if ~P.instantSludgeRemoval
        dVdt = P.dVdt_pu;
        if P.variableInGrams 
            dXdt = dXdt+P.dX_g_pu;
        else
            dX_c_pu = (P.dX_g_pu*V-P.dVdt_pu*X_g)/(V^2); % cdot=(dot(m/V)=dotm*V-dotV*m)/V*V % this is not global because V and X_g change at every round
            dXdt = dXdt+dX_c_pu;
            dS_g_pu = zeros(size(S_g)); % Subsgrams does not change during the purge, there is no efflux of solubles
            dSdt = dSdt+(dS_g_pu*V-dVdt*S_g)/(V^2); % cdot=(dot(m/V)=dotm*V-dotV*m)/V*V
        end
    else
        % do nothing, this is already done in defineInitialConditions_intermediateSteps
    end
 end   
end

if P.onEfflux && P.influxActive % here NO biomass is removed. influxActive because I make sure to take out only if I put in
 dVdt = P.dVdt_ef; %QEffl = Q2
 if P.variableInGrams % only Subsgrams change, they are removed by the efflux
    dSdt = dSdt+dVdt*S_c;
 else % concentration of Subs does not change, removed with the water, only the volume changes
    dX_g_ef = zeros(size(dXdt));
    dXdt = dXdt+(dX_g_ef*V-dVdt*X_g)/(V^2); % cdot=(dot(m/V)=dotm*V-dotV*m)/V*V
 end
end
 
if P.onAer % in this case I have that the initial concentration of oxygen equal to the targetO2conc and the derivatives remains zero
 dSdt(P.sO2) = 0;
end

dxdt(1:P.nPart) = dXdt;
dxdt((P.nPart+1):(P.nPart+P.nSolub)) = dSdt; % adding differentials for nutrients
dxdt(end) = dVdt;
return


function P = checksInterm(P)
% Intermediate checks during the ODEs integration

S = P.this_x_all_0((P.nPart+1):(P.nPart+P.nSolub));
V = P.this_x_all_0(end); 

% - Check oxygen concentration
MEstr = 'check how oxygen is set during the cycles';
if P.onAer && (((~P.variableInGrams && ~reallyEqual(S(P.sO2),P.O2_c_ref))) || (P.variableInGrams && ~reallyEqual(S(P.sO2)/V,P.O2_c_ref))) 
 ME = MException(MEstr);
 throw(ME)
end

% - Check that the Volume changes are correct at the end of the phase. 
% NB there might be differencies due to of numerical inaccuracies (beacause QEfflux and QPurge are calculated based on X(t))
switch P.pastPhase
 case 'onFeed'
    indPhace1Cy = mod(P.indPhase,P.nPhases_1C);
    nFeed_1Cy = sum(P.phaseSeq_1C==3); % I have 3 feed phases per cycle
    if ~reallyEqual(V,(P.V_0+P.V_in_1Cy/nFeed_1Cy*(indPhace1Cy-1)/nFeed_1Cy)) 
      ctrl = P.V_0+P.V_in_1Cy/nFeed_1Cy*(indPhace1Cy-1)/nFeed_1Cy;
      P.epsyFeed = max(P.epsyFeed,abs(V-ctrl));
      V = ctrl;
    end
 case 'onPurge'
    if P.washOutParticles
        if ~reallyEqual(V,(P.V_0+P.V_in_1Cy-P.V_pu_NextCy))
          ctrl = (P.V_0+P.V_in_1Cy-P.V_pu_NextCy);
          P.epsyPurge = max(P.epsyPurge,abs(V-ctrl));
          V = ctrl;
        end
    end
 case 'onEfflux'
    if P.influxActive
        if P.indPhase>1
            MEstr = 'at the end of one cycle the reactor volume should be brought back to its equilibrium value';
            reallyEqual(P.V_0,(P.V_0+P.V_in_1Cy-P.V_pu_NextCy-P.V_ef_NextCy),MEstr);
            [~,V,~] = reallyEqual(V,P.V_0,'~');
        end
    end
end
if P.phaseCtrl == 1 && ~reallyEqual(V,P.V_0)
 P.exche = [P.exche; V];
 ME = MException('at the end of one cycle the reactor volume should be brought back to its equilibrium value');
 throw(ME)
end
return

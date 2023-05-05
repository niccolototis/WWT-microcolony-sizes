function [grNet] = compute_mu(S_c,sp,P) 
 switch sp
    case 'HET'
        % I am not using it anyway, in this study the PBE models only
        % AUT microcolonies
        grNet = P.MuH.*(S_c(:,P.sS)./(P.Ks+S_c(:,P.sS))).*(S_c(:,P.sO2)./(P.KOH+S_c(:,P.sO2))).*(S_c(:,P.sNH4)./(P.KNH_mod+S_c(:,P.sNH4)))+ ... % aerobic growth
                P.MuH.*P.Etag.*(S_c(:,P.sS)./(P.Ks+S_c(:,P.sS))).*(P.KOH./(P.KOH+S_c(:,P.sO2))).*(S_c(:,P.sNO)./(P.KNO+S_c(:,P.sNO))).*(S_c(:,P.sNH4)./(P.KNH_mod+S_c(:,P.sNH4)))+ ... % anoxic growth
                -P.bH;
            
    case 'AUT'
        grNet = P.MuA.*(S_c(:,P.sNH4)./(P.KNH+S_c(:,P.sNH4))).*(S_c(:,P.sO2)./(P.KOH+S_c(:,P.sO2)))+...
        -P.bA;
 end
return

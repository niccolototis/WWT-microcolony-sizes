function dSdt_aSp = updateSubs(Subs,Vg_Pos,Vg_Neg,aSpN,x2_widths)
 global HET AOB NOB  NO2 NO3 S O2 NH4 decay EPS modelEPS I decayActive Y S_in modelI pres
 dSdt_aSp = zeros(size(S_in));
 switch aSpN
    case 'HET'
        dEPS = 0;
        if decayActive && modelEPS  % remember EPS only produced ny HET
            dEPS = Subs(EPS)*decay(EPS);
            dSdt_aSp(EPS) = Y(EPS)/Y(HET)* sum(Vg_Pos.HET.allModes.*x2_widths)*pres(HET)-dEPS;
        end
        if modelI 
            dSdt_aSp(I) = +Y(I) * sum(Vg_Neg.HET.*x2_widths)*pres(HET);
        end
        dSdt_aSp(S) = -1/Y(HET) * sum(Vg_Pos.HET.allModes.*x2_widths)*pres(HET)+dEPS;
        dSdt_aSp(O2) = -(1-Y(HET)-Y(EPS))/Y(HET) * sum(Vg_Pos.HET.AER.*x2_widths)*pres(HET);
        dSdt_aSp(NO2) = -(1-Y(HET)-Y(EPS))/(1.71*Y(HET)) * sum(Vg_Pos.HET.NO2.*x2_widths)*pres(HET);
        dSdt_aSp(NO3) = -(1-Y(HET)-Y(EPS))/(2.86*Y(HET))* sum(Vg_Pos.HET.NO3.*x2_widths)*pres(HET); %% here I see that NO3 when is metabolized by HET.NO3 either goes to HET.NO3 biomass, either to            
    case 'AOB'
        if modelI 
            dSdt_aSp(I) = +Y(I) * sum(Vg_Neg.AOB.*x2_widths)*pres(AOB);
        end
        dSdt_aSp(S) = +(1-Y(I)) * sum(Vg_Neg.AOB.*x2_widths)*pres(AOB);
        dSdt_aSp(O2) = -(3.42-Y(AOB))/Y(AOB) * sum(Vg_Pos.AOB.allModes.*x2_widths)*pres(AOB);
        dSdt_aSp(NH4) = -1/Y(AOB)*sum(Vg_Pos.AOB.allModes.*x2_widths)*pres(AOB);
        dSdt_aSp(NO2) = +1/Y(AOB) * sum(Vg_Pos.AOB.allModes.*x2_widths)*pres(AOB);
    case 'NOB'
        if modelI 
            dSdt_aSp(I) = +Y(I) * sum(Vg_Neg.NOB.*x2_widths)*pres(NOB);
        end
        dSdt_aSp(S) = +(1-Y(I)) * sum(Vg_Neg.NOB.*x2_widths)*pres(NOB);
        dSdt_aSp(O2) = -(1.15-Y(NOB))/Y(NOB) * sum(Vg_Pos.NOB.allModes.*x2_widths)*pres(NOB);
        dSdt_aSp(NO2) = -1/Y(NOB)*sum(Vg_Pos.NOB.allModes.*x2_widths)*pres(NOB);
        dSdt_aSp(NO3) = +1/Y(NOB)*sum(Vg_Pos.NOB.allModes.*x2_widths)*pres(NOB); 
 end
return

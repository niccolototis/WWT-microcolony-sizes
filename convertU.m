function [] = convertU()
global solPartic solubles q_O2 S_star_O2  diamMin_T diamInitPart_T diamMax_T 
global mu_max K S_in decay      newOneMeter   HRT SDT   grOneKg miColDensity 
 %conversions
 miColDensity = miColDensity*(1/grOneKg)*(1e03)*(1e03); % from g/mL to Kg/m3:
 diamMin_T = diamMin_T  *newOneMeter;
 diamInitPart_T = diamInitPart_T  *newOneMeter;
 diamMax_T = diamMax_T  *newOneMeter; 
 % ------------
 q_O2 = q_O2*(grOneKg/(hoursOneDay*newOneMeter^3));  %g/(h*micron^3)      
 S_star_O2 = S_star_O2*(grOneKg/newOneMeter^3); %g/micron^3      
 HRT = HRT*hoursOneDay; %h                      
 SDT = HRT*(alpha+beta)/(beta*(1+alpha)); %h 
 S_in(solubles) = S_in(solubles) .*(grOneKg/newOneMeter^3); %g/micron^3 @@NT
 S_in(solPartic) = S_in(solPartic).*grOneKg; %@@NT . NB different units in S_in!!
 K = K .*(grOneKg/newOneMeter^3); %g/micron^3 @@NT
 decay = decay .*(1/hoursOneDay); %1/h
 mu_max = mu_max .*(1/hoursOneDay); %1/h
return

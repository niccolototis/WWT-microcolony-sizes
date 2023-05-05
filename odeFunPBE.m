function dN_all_aSp_aChar = odeFunPBE(t,N_all_aSp_aChar,aStp,BrRateKer_B,BrRateKer_V,t_all_ASM1seq,mu_ASM1seq_tall,aChar,sp,P)  
global steP step
% In this code it is assumed that the whole bag model has already been run
% and that the trajectories for particles solubles and vol have already been
% computed

if ~P.identify_lambda && ~P.sweep_lambda
steP=steP+1; step=step+1; 
end

% disp(t)
N_all_aSp_aChar= fixLowerThreshold_fast(N_all_aSp_aChar,P); 

% Take the solubles at the tpoint lower or equal to t (not go forward in time)
t_diffs=t-t_all_ASM1seq;
t_diffs(t_diffs<0)=Inf;
[~,idxMin]=min(t_diffs);
P.thisTpoint=idxMin;

mu_ASM1seq=mu_ASM1seq_tall(P.thisTpoint);

%%%%%%%%%%%%%%%####################################%%%%%%%%%%%#############

Nti_aChar=N_all_aSp_aChar(1:P.nCells);
x1=N_all_aSp_aChar(P.nCells+1); 
w1=N_all_aSp_aChar(P.nCells+2); 
[dx1,dw1]=deal(0);
if length(N_all_aSp_aChar)~=(P.nCells+2)
    ME = MException('check dimensions of arrays');
    throw(ME)
end

% if  mod(steP,1000)==1 
% %     M_NisTilde=true;
%     Marg_x1_aChar=computeMarg(Nti_aChar,P.M_NisTilde,w1,aStp.x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
%     disp(['t: ',num2str(t),' nChar:',num2str(aChar),'  aSp: ',char(P.spNames{sp}),'  Marg_x1_thisChar: ', num2str(Marg_x1_aChar)]);
% end


% +++++++++++++++++ MOC: CHARACTERISTIC ODES +++++++++++++++++++
% Characteristics ODEs used to compute dx1 dw1
%     growth_x1=; %this is div(f(xi(t))) in eq 4-5 of paper. NB should not depend on x2! (it's an additional assumption besides being non growth-related)
%     dx1 = x1*growth_x1;
%     dw1 = w1*growth_x1;


% +++++++++++++++++ START CAT ++++++++++++++++++++++++++++++++++
[B,V,D]=deal(zeros(size(Nti_aChar)));

[Bin,Vin] = FlowIn(Nti_aChar,P);
checkNNG(Bin);
checkNNG(Vin);
B=B+Bin;
V=V+Vin;

if P.implementGrowth
    [Bg,Vg,Dg] = GrowthPBE(Nti_aChar,aStp,mu_ASM1seq); 
    checkNNG(Bg);
    checkNNG(Vg);
    B=B+Bg;
    V=V+Vg;
    D=D+Dg;
end

if P.implementBreakage
    Nti_aChar_toBreak=skip1stRef(Nti_aChar,P.FirstRefNoBreak); 
    [Bb,Vb,Db] = BreakagePBE_WWT(Nti_aChar_toBreak,@(x2_all,aSp) BrRateWWT(aStp.x2_all,P),BrRateKer_B,BrRateKer_V,aStp,P); % BreakagePBE(N,BrRate,BrRateKer_B,BrRateKer_V)
    checkEqualStr_ForLoop(sum(Bb),2*sum(Db),'breakage',P);
    checkNNG(Bb);
    checkNNG(Db);
    B=B+Bb;
    D=D+Db;
    V=V+Vb;
end

[Dw] = washOut(Nti_aChar,P);
checkNNG(Dw);
D=D+Dw;
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

checkNumericalInaccuracy(B,D,'PRE',P)
B=setNNG(B);
V=setNNG(V);
B_CA = computeRearrangeCAT_new(B,V,aStp,P,t); %  computeRearrangeCAT(B,V,x2_all,x2_lb,x2_ub,t);
checkNumericalInaccuracy(B_CA,D,'POST',P)   
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dNti_aChar=B_CA-D; % Death does not affect the redistribution among different cells in the CAT        
dN_all_aSp_aChar=[dNti_aChar;dx1;dw1];
% disp(strcat(['Minimum dNdt>0 : ', num2str(min(dN_all_aSp_aChar(dN_all_aSp_aChar>0))), ' - Max dNdt>0 : ', num2str(max(dN_all_aSp_aChar(dN_all_aSp_aChar>0))) ]))
return

function [Bin,Vin] = FlowIn(Nti_aChar,P)
[Bin,Vin]=deal(zeros(size(Nti_aChar)));
if P.particleInflux && P.onFeed %onPurge 
    %define 
%     Bin=N*(1/SRT)*(1/cyclesOneDay); 
end
return


function Dw = washOut(N,P)
Dw=0;
switch P.modelWaterSludgeSegreg 
case 'segregated'
    if P.onPurge && P.washOutParticles
        if ~P.instantSludgeRemoval 
            Dw=P.Dw_purge;
        else
            %if instantSludgeRemoval do nothing, this is already done in defineInitialConditions_intermediateSteps
        end    
    end
case 'stirred'
    %             if onEfflux % here also the biomass is removed (paseCtrl=7). I do not have a separate purge phase
    %             dVdt=-QEffl; %QEffl=Q2
    %             if variableInGrams %NB Subs are conc, X are grams
    %                 if ~instantSludgeRemoval
    %                     dXdt=dXdt-QEffl*(X/V); % cdot=(dot(m/V)=dotm*V-dotV*m)/V*V
    %                 else
    %                     % do nothing, this is already done in defineInitialConditions_intermediateSteps
    %                 end
    %                 dSdt=dSdt-QEffl*Subs;
    %             else
    %                 if instantSludgeRemoval
    %                     mX=X*V; %compute the number of grams of particulate components
    %                     dXdt=dXdt-(mX*dVdt)/(V*V); % cdot=(dot(m/V)=dotm*V-dotV*m)/V*V
    %                 end
    %             end
    %         end
end

return

function []=checkNumericalInaccuracy(B,D,PREoPOST,P)
global epsylon
    if ~P.identify_lambda && ~P.sweep_lambda
        if ~P.implementGrowth && P.implementBreakage && ~P.washOutParticles && ~P.influxActive
            tepsy=abs(sum(B)-2*sum(D));
            if tepsy>epsylon % they should sum up to one aprt from where I have no births. For some numerical reason there is a little inaccuracy, so a epsylos should be considered/ here epsylon=1e-15
                epsylon=tepsy;
                disp(['epsylon--',PREoPOST, num2str(tepsy)])
                %     ME = MException('total B and total D should sum up to the same number');
                %     throw(ME)
            end
        end    
    end
return

function [Nti_aChar]= skip1stRef(Nti_aChar,FirstRefNoBreak)
% This function puts the number density of the first reference to zero if
% the first reference coincides with the size of the smallest particle that
% can be create through breakage/erosion (== x2SCell in this case).
% In this way particles with size x2=x2SCell do not erode -they would
% erode into themselves, but when nu>1 this creates problems-
% NB: this has been already taken into account also in the construction of
% BrRateKer_B and BrRateKer_V
if FirstRefNoBreak
    Nti_aChar(1)=0;
end
return


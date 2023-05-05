function dN_all = odeFun_2D_PBE(t,N_all,x2_all,x2_b,x2_lb,x2_ub,x2_x0,N_x0,x2_widths,BrRateKer_b,BrRateKer_v,ghost1,ghostEnd)
global nCells nChar nCombo steP implementBreakage implementGrowth thr p aChar M_NisTilde

steP = steP+1;
N_all = fixLowerThreshold(N_all,thr);
Nti_aChar = N_all(1:nCells);
Nti_aChar = ghostToZero(Nti_aChar);
x1_aChar = N_all(nCells+1);
w1_aChar = N_all(end);
dx1 = zeros(size(x1_aChar));
dw1 = zeros(size(w1_aChar));
dNti = zeros(size(Nti_aChar));

if mod(steP,50) ==1 
% M_NisTilde = true;
 Marg_x1_aChar = computeMarg(Nti_aChar,M_NisTilde,w1_aChar,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
 disp(['t: ',num2str(t),' nChar:',num2str(aChar),' maxN: ',num2str(max(Nti_aChar)),'  Marg_x1_thisChar: ', num2str(Marg_x1_aChar)]);
end
% +++++++++++++++++ MOC: CHARACTERISTIC ODES +++++++++++++++++++
% Characteristics ODEs used to compute dx1_all dw1_all
% !!NB should NOT depend on x2 coord, otherwise using MOC does not
% simplify the discretization wrt a CAT_2D!
% dx1 = x1*(p.a-p.b);%*x2_all 
% dw1 = w1*(p.a-p.b);%*x2_all);

% +++++++++++++++++ START CAT ++++++++++++++++++++++++++++++++++
[B,V,D] = deal(0);
if implementGrowth
 [Bg,Vg,Dg] = Growth(Nti_aChar,x2_all,N_x0,x2_x0,x1_aChar);
 checkEqualStr(Bg,Dg,'growth'); 
 checkNNG(Bg);
 checkNNG(Dg);
 checkNNG(Vg);
 B = B+Bg;
 D = D+Dg;
 V = V+Vg;
end
if implementBreakage
 [Bb,Vb,Db] = BreakagePBE(Nti_aChar,@(x2_all,aSp) BrRate(x2_all,[]),BrRateKer_b,BrRateKer_v,x2_all,x2_widths,[]);
 checkEqualStr(sum(Bb),2*sum(Db),'breakage');
 checkNNG(Bb);
 checkNNG(Db);
 B = B+Bb;
 V = V+Vb;
 D = D+Db;
end
checkNumericalInaccuracy(B,D,'PRE')
if ~(isequal(B,0)&&isequal(V,0)&&isequal(D,0))
 B_CA = computeRearrangeCAT(B,V,x2_all,x2_lb,x2_ub,t,ghost1,ghostEnd);
%     checkNumericalInaccuracy(B_CA,D_CA,'POST')
 B_CA = ghostToZero(B_CA);
 D = ghostToZero(D);
 dNti = B_CA-D; 
 dNti = ghostToZero(dNti);
end
dN_all = [dNti;dx1;dw1];
return

function [Bg,Vg,Dg] = Growth(N,x2_all,N_x0,x2_x0,x1)
global steP
 gr = grRate(x2_all,x1);
 [grPos,grNeg] = deal(zeros(size(gr))); % need to split apart the positive vs negative growth
 grPos(gr>= 0)=gr(gr>=0);
 grNeg(gr<0) = abs(gr(gr<0));
 Bg_noDecay = N * N_x0 .* grPos; %?Aerobic growth of heterotrophs HET
 Bg_decayRate = N * N_x0 .* grNeg;
 Bg = Bg_noDecay  +   Bg_decayRate; %impo the plus sign here!!! sono comunque nascite, che coincidono con le morti. IL DECAY O DEGROWTH EMERGE DAL FATTO CHE LE NEW PARTICLES GENERATE BY DECAY HANNO VOLUME (v+x2_x0)!
 Dg = Bg;
 Vg = Bg_noDecay.*(x2_all+x2_x0)  +   Bg_decayRate.* (x2_all-x2_x0); 
 % I wanna have zero growth in ghost cells
 Bg = ghostToZero(Bg);
 Dg = ghostToZero(Dg);
 Vg = ghostToZero(Vg); 
return



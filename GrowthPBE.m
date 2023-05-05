function [Bg,Vg,Dg] = GrowthPBE(N,aStp,mu_ASM1seq)
grNet = mu_ASM1seq*aStp.x2_all; % xdot=mu*x
% --
[grPos,grNeg] = deal(zeros(size(grNet))); % need to split apart the positive vs negative growth. @NB as mu_ASM1seq is size-independent, the entries of the vector grNet are all neg or all pos
grPos(grNet>= 0)=grNet(grNet>=0);
grNeg(grNet<0) = abs(grNet(grNet<0));

% --
Bg_Pos = N * aStp.N_x0 .* grPos; %?Aerobic growth of heterotrophs HET
Bg_Neg = N * aStp.N_x0 .* grNeg;
Bg = Bg_Pos + Bg_Neg; %impo the plus sign here!!! sono comunque nascite, che coincidono con le morti. IL DECAY O DEGROWTH EMERGE DAL FATTO CHE LE NEW PARTICLES GENERATE BY DECAY HANNO VOLUME (v+x2_x0)!
Dg = Bg;
Vg_Pos = Bg_Pos .* (aStp.x2_all+aStp.x2_x0);
Vg_Neg = Bg_Neg .* (aStp.x2_all-aStp.x2_x0);
Vg = Vg_Pos + Vg_Neg ; 

%     Exception of SCell decaying ==> I am solving this issue in
%     computeCAT_new function
% [Bg,Vg,Dg] = SCellDecay(aStp,Bg_Neg,Vg_Neg,Bg,Vg,Dg);

% I wanna have zero growth in ghost cells
% [Bg(ghostIdx),Dg(ghostIdx),Vg(ghostIdx)] = deal(0);
return

% function [Bg,Vg,Dg] = SCellDecay(aStp,Bg_Neg,Vg_Neg,Bg,Vg,Dg)
% % This function considers the exception of the first cell grid: SCells
% % cannot shrink below their size, thus if grNet is negative, particles are
% % removed instead of shrinked. The number of removed particles (Dg) is computed
% % from the total particle volume to be shrinked
% global noGhostInd
% ind1Ref = noGhostInd(1);
% if Vg_Neg(ind1Ref)>0
%     VgLost1Ref = Bg_Neg(ind1Ref)*aStp.x2_all(ind1Ref)-Vg_Neg(ind1Ref); % from equations 42 and 40
%     Bg(ind1Ref) = 0;
%     Vg(ind1Ref) = 0;
%     Dg(ind1Ref) = VgLost1Ref/aStp.x2_all(ind1Ref);
% end
% return

% This is the reference code used for 2D_PBE used foc CDC paper
% global steP
% gr = grRate(x2_all,x1);
% [grPos,grNeg] = deal(zeros(size(gr))); % need to split apart the positive vs negative growth
% grPos(gr>= 0)=gr(gr>=0);
% grNeg(gr<0) = abs(gr(gr<0));
% Bg_noDecay = N * N_x0 .* grPos; %?Aerobic growth of heterotrophs HET
% Bg_decayRate = N * N_x0 .* grNeg;
% Bg = Bg_noDecay   +   Bg_decayRate; %impo the plus sign here!!! sono comunque nascite, che coincidono con le morti. IL DECAY O DEGROWTH EMERGE DAL FATTO CHE LE NEW PARTICLES GENERATE BY DECAY HANNO VOLUME (v+x2_x0)!
% Dg = Bg;
% Vg = Bg_noDecay.*(x2_all+x2_x0)   +   Bg_decayRate.* (x2_all-x2_x0); 
% % I wanna have zero growth in ghost cells
% Bg = ghostToZero(Bg);
% Dg = ghostToZero(Dg);
% Vg = ghostToZero(Vg); 
% return

function [Bb,Vb,Db] = BreakagePBE_WWT(N,BrRate,Ker_B,Ker_V,aStp,P)
% BrRate has to be a function handle
% Bb = nx1; Ker_B=nxn; (x2_widths .* BrRate(v,aSp) .* N)=nx1; "Ker_B *" is a
% matrix product, I am summing up along the rows of Ker_B
% % EITHER THIS VERSION IS CORRECT
%     Db = BrRate(aStp.x2_all,aSpPBE) .* N;     
%     % I multiply each column of by the corrisponding element of Ker_B
%     % and Ker_V by the vector (x2_widths .* BrRate(v,aSp) .* N)
%     Bb = nu * Ker_B * (aStp.x2_widths .* BrRate(aStp.x2_all,aSpPBE) .* N); 
%     Vb = nu * Ker_V * (aStp.x2_widths .* BrRate(aStp.x2_all,aSpPBE) .* N);
%     Db = ghostToZero(Db); 
%     Bb = ghostToZero(Bb); 
%     Vb = ghostToZero(Vb); 
%     
% OR THE PREVIOUS VERSION
    Db = BrRateWWT(aStp.x2_all,P) .* N;     % BrRate is a function handle
    Bb = P.nu * Ker_B * (BrRate(aStp.x2_all) .* N); % NB here there should NOT be the factor 2BrRate(2*v,aSp), because nnzero entries of BrRateKer_b already select the specific N(mother_idx) and BrRate(v(mother_idx),aSp) where v(mother_idx) are those whose breakage fall in the vD interval considered
    Vb = P.nu * Ker_V * (BrRate(aStp.x2_all) .* N); 
return


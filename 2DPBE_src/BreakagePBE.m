function [Bb,Vb,Db] = BreakagePBE(N,BrRate,BrRateKer_b,BrRateKer_v,v,v_widths)
global nu
%     Db = BrRate(v,aSp) .* N;     % BrRate is a function handle
%     Bb = nu * BrRateKer_b * (BrRate(v,aSp) .* N); % NB here there should NOT be the factor 2BrRate(2*v,aSp), because nnzero entries of BrRateKer_b already select the specific N(mother_idx) and BrRate(v(mother_idx),aSp) where v(mother_idx) are those whose breakage fall in the vD interval considered
%     Vb = nu * BrRateKer_v * (BrRate(v,aSp) .* N); 
%     Db = ghostToZero(Db); %OPTIONAL
%     Bb = ghostToZero(Bb); %OPTIONAL
%     Vb = ghostToZero(Vb); %OPTIONAL
    
    
    
    % EITHER THIS VERSION IS CORRECT
    Db = BrRate(v,aSp) .* N;     
    % I multiply each column of by the corrisponding element of Ker_B and Ker_V by the vector (x2_widths .* BrRate(v,aSp) .* N)
    Bb = nu * BrRateKer_v * (v_widths .* BrRate(v,aSp) .* N); 
    Vb = nu * Ker_V * (v_widths .* BrRate(v,aSp) .* N); 
%     Bb = nu * Ker_B * ( BrRate(v,aSp) .* N); 
%     Vb = nu * Ker_V * ( BrRate(v,aSp) .* N); 
    Db = ghostToZero(Db); %OPTIONAL
    Bb = ghostToZero(Bb); %OPTIONAL
    Vb = ghostToZero(Vb); %OPTIONAL
return


function [moment] = computeMoments(M_N_0,M_NisTilde,x2_all,x2_widths,x1_all,w1_all,nChar,momentForWhatVariable,numberOfMoment)
% Ref: 
% [1] ?STEFFEN WALDHERR, PHILIP TRENNT, AND MUBASHIR HUSSAIN?Hybrid Simulation of Heterogeneous Cell Populations
% [2] http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Importance_Sampling.html

% Calculates the moments as described in [1]
 if isequal(numberOfMoment,0)
    x = ones(size(x1_all)); % ininfluent
 else
    switch momentForWhatVariable
        case 'x1'
            x = x1_all;
        case 'x2'
            x = x2_all;
    end
 end
 M_N_0_x = M_N_0.*(x.^numberOfMoment);
 Marg_x1 = computeMarg(M_N_0_x,M_NisTilde,w1_all,x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
 Marg_x1_ti = Marg_x1.*w1_all; % Fundamental, need to return to tilde before the next line
 moment = computeMarg(Marg_x1_ti,M_NisTilde,[],[],nChar,'x2');
return

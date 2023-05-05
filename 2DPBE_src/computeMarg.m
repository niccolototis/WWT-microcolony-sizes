function [Marg] = computeMarg(M_N_0,M_NisTilde,w1_all,x2_widths,nx1,MargForWhat)
% Calculates the Marginal over x1 or x2 from a joint produced with parameters mu and sigma
% M_N_0 = computeJoint(x1pts,x2pts,nx1,nx2, mu, sigma);
% @@NB; This function returns always a marginal which is weight-independent! Mind this when applying this
% function twice consequently

% Ref: 
% [1] ?STEFFEN WALDHERR, PHILIP TRENNT, AND MUBASHIR HUSSAIN?Hybrid Simulation of Heterogeneous Cell Populations
% [2] http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Importance_Sampling.html

switch MargForWhat
 case 'x1' 
    sumAcross = 1; % sums vertically (along the rows) all the terms in each column
    % The marginal for x1 is computed taking the integral over the x2 feature space. This integral is computed in the classical discretized form sum(fx.*dx)
    Marg = sum((x2_widths.*M_N_0),sumAcross); % 
    if M_NisTilde
        Marg = Marg./w1_all; % see Table 1 in [1]
    end
 case 'x2' % The marginal for x2 is computed taking the integral over the x1 feature space. This integral is computed with MC technique
    if ~M_NisTilde
        ME = MException('I can only compute the marginal for x2 (integrating over x1 coordinate) starting from n_tilde as in [1], because I do not have x1_widths');
        throw(ME) 
    end
    M_Nti_0 = M_N_0;
    sumAcross = 2;% sums horizontally (along the columns) all the terms in each row
    Marg = 1/nx1*sum(M_Nti_0,sumAcross); % here I do not know the widths, but I am performing the importance sampling. See "Marginal for the y variables" in Table 1 in [1]
end
return

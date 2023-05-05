% function [x2_all,x2_b,x2_lb,x2_ub,x2_x0,N_x0,x2_widths,ghost1,ghostEnd, M_N_0] = defineIC_SizeGRate_forLoop()
function [M_Nti_0_full,x2_all,x2_b,x2_lb,x2_ub,x2_x0,N_x0,x2_widths,ghost1,ghostEnd] = defineIC_2D_PBE()
global nCells nChar steP noGhostInd implementGrowth implementBreakage thr p decayRateActive epsylon lostFirstMoment ghostIdx nCombo plot3Darmin DivisionRate
global twoPointsFormula mu_x1 mu_x2 sigma_x1 sigma_x2 x2Min x2Max x1Min x1Max nD x1Sampling nCellsGh Ntot_0

% Ref: 
% [1] ?STEFFEN WALDHERR, PHILIP TRENNT, AND MUBASHIR HUSSAIN?Hybrid Simulation of Heterogeneous Cell Populations
% [2] http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Importance_Sampling.html

% OUTPUTS: 
% M_Nti_0_full is the matrix containing as many rows as nCells and as many columns as nChar. Each column
% of this matrix provides a initial condition for the PBE modeling the evolution of Ñ in time, as
% described in [1]

% scaleDots = 'g*';
twoPF = {'false2PF'};
if twoPointsFormula && implementGrowth
 twoPF = ['true2PF',twoPF];
else
 twoPointsFormula = false;
end

%%%%%%%%%%%%%%%%%%%%%%% discretization attribute space %%%%%%%%%%%%%%%%
% figure

mu = [mu_x1 mu_x2]; % in cm
sigma = [sigma_x1 sigma_x2].*eye(2); % in cm


% 1 - discretized grid of growth related coordinate x2
x2_b = linspace(x2Min,x2Max,nCells+1)';
x2_lb = x2_b(1:end-1);
x2_ub = x2_b(2:end);
x2_all = (x2_lb+x2_ub)/2; % reference at the midpoint
[x2_b,x2_lb,x2_ub,x2_all,x2_widths,ghost1,ghostEnd] = addGhostCells(x2_b,x2_all);
nCombo = nChar*nCells;

m_0 = 1; % complexive volume of cells of minimal size that are used to compute the growth term 
x2_x0 = min(x2_all(x2_all>0))/100; % used to implement growth with the CAT. size of the "diffused particles"
N_x0 = m_0/x2_x0; % used to implement growth with the CAT. ndf of the "diffused particles"


% 2: define coordinates for x1
figure;
switch x1Sampling
 case 'MH-sampling'
    % 2.a - use Metrppolis Hastings algorithm for MC sampling from the marginal distribution (of the joint) for non-growth related coordinate x1 (integrating over all x2_all points)
    % smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd)
    start = mu_x1;                    
    discarded = 200;
    nsamples = nChar;
    nx1 = 1;
    MargForWhat = 'x1';
    thisPdf = @(x1pt) computeMarg(computeJoint(x1pt,x2_all,nx1,nCellsGh,mu,sigma),[],x2_widths,nx1,MargForWhat); %
    delta = sigma_x1;
    proppdf = @(x,y) normpdf(x,y,delta); %The proposal distribution q(x,y) gives the probability density for choosing x as the next point when y is the current point. Use here current point y as the mean of the proposal distrib function == normpdf
    proprnd = @(y) y + randn*delta;  % proprnd defines the random number generator for the proposal distribution (qui genero x a seconda di y
    [x1_all,accept] = mhsample(start,nsamples,'pdf',thisPdf,'proppdf',proppdf, 'proprnd',proprnd,'symmetric',false,'burnin',discarded);
    x1_all = sort(x1_all);
    h = histfit(x1_all,50);
    hold on
    plot(x1_all,normpdf(x1_all,mu_x1,sigma_x1),'LineWidth',2);
    h(1).FaceColor = [.8 .8 1];
 case 'normal'
    % 2.b -  sampling from the Normal distribution for non-growth related coordinate x1
    x1_all = normrnd(mu_x1,sigma_x1,nChar,1);
 case 'uniform'
    % 2.c -  linear sampling for non-growth related coordinate x1
    x1_all = linspace(x1Min,x1Max,nChar);
end
x1_all = sort(x1_all);


% 3 - calculate the initial marginal frequency for the chosen x1_all, calculate associated weights 
M_N_0 = computeJoint(x1_all,x2_all,nChar,nCellsGh, mu, sigma); %IMPO! this x1_all needs to be here
M_N_0 = M_N_0* Ntot_0; % !! IMPO @@NT

MargForWhat = 'x1'; % function [Marg]=computeMarg(M_N_0,w1_all,x2_widths,nx1,MargForWhat)
M_NisTilde = false;
Mx1_0 = computeMarg(M_N_0,M_NisTilde,[],x2_widths,nChar,MargForWhat);
%         MC integration of marginal Mx1_0 http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Importance_Sampling.html
int_Joint_0 = 1* Ntot_0; % e This is the integral of the joint == integral of the marginal for x1 over the same coordinate x1 .
w1_all = (int_Joint_0./Mx1_0);
px1_all = (1./w1_all);
if ~isequal(num2str((1/nChar)*sum(Mx1_0.*w1_all)),num2str(int_Joint_0))
 ME = MExeption('something is going wrong with the importance sampling, look to the material here below');
 throw(ME);
end
%{
Following [2]
I = integral to be estimated with Monte Carlo method
I = 1/N*symsum( f(x(i))/p(x(i)) ),i,0,N)
p(x(i)) = c*f(x(i)) % I want the probability of my samples to have more or less the same distribution of
the pdf I want to calculate the integral of
by normalization c = 1/int(f(x)0
Since c is a constant, each estimate has the same value, and the variance is zero! Of course, this is ludicrous since we wouldn?t bother using Monte Carlo if we could integrate f(x) directly. However, if a density p(x) can be found that is similar in shape to f(x), variance decreases.
%}

% 4 - calculate Ñ as in equation 10 in [1]
M_Nti_0 = M_N_0.*w1_all; % Each column multiplied by a specific weight
MargForWhat = 'x2';
M_NisTilde = true;
Mx2_0 = computeMarg(M_Nti_0,M_NisTilde,[],x2_widths,nChar,MargForWhat);

% 5 - M_N_0 is the joint probability. 
% compare the 2D integrals calculated inn the two following ways. Fon nChar -> inf they should coincide
M_NisTilde = true;
Marg_x1 = computeMarg(M_Nti_0,M_NisTilde,w1_all,x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
Marg_x1_ti = Marg_x1.*w1_all; % Fundamental, need to return to tilde before the next line
check2DIntegral_v1 = computeMarg(Marg_x1_ti,M_NisTilde,[],[],nChar,'x2');
% vs
Marg_x2 = computeMarg(M_Nti_0,M_NisTilde,[],[],nChar,'x2');
M_NisTilde = false;
check2DIntegral_v2 = computeMarg(Marg_x2,M_NisTilde,[],x2_widths,[],'x1');
% 
if ~isequal(num2str(check2DIntegral_v1),num2str(check2DIntegral_v1)) || ~isequal(num2str(check2DIntegral_v1),num2str(int_Joint_0))
 ME = MException('These integrals should match');
 throw(ME) 
end 
if (abs((check2DIntegral_v1-Ntot_0)/Ntot_0)> 1e-6) || (abs((check2DIntegral_v1-Ntot_0)/Ntot_0)> 1e-6)
 disp(['v1:', num2str(check2DIntegral_v1),' - v2:', num2str(check2DIntegral_v2), ' Ntot_0:',num2str(Ntot_0)])
 ME = MException('I want the integral of the jont ndf to equal the specified initial number of particles');
 throw(ME) 
end       
% N = M_Nti_0(:); % this are all initial Ñ for all charachteristics
M_Nti_0_full = [M_Nti_0; x1_all; w1_all]; 

%%%%%%%%%%%%%%%%%%%%%%% plot IC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(x1_all, Mx1_0, 'LineWidth', 2)
title('Marginal for x1')
% figure;
% plot(-Mx1_0,x1_all, 'LineWidth', 2)
% hold on
% plot(range(Mx1_0), x1_all, 'LineWidth', 2,'Color','r')
figure;
plot(x2_all(noGhostInd), Mx2_0(noGhostInd), 'LineWidth', 2) % remove ghost cells
title('Marginal for x2')
epsylon = 0; % This is used to evaluate numerical inaccuracy
steP = 0;
return






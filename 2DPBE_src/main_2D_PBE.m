function [] = main_2D_PBE()
clear all
close all

global nCells nChar step noGhostInd implementGrowth implementBreakage thr p epsylon lostFirstMoment ghostIdx nCombo plot3Darmin BreakageRate
global k twoPointsFormula mu_x1 mu_x2 sigma_x1 sigma_x2 x2Min x2Max x1Min x1Max aChar g nCellsGh Ntot_0 M_NisTilde approx_lb_max approx_ub_max
global epsyBreak epsyGrow nInconsGrowth nInconsBreak lostFirstMoment_lb lostFirstMoment_ub tFlaglb tFlagub nD x1Sampling GrowthRate

rmpath('/Users/niccolo/git_KU/benelux2019/CAT/CAT_NT/CAT_steffenspostdoc')
rmpath('/Users/niccolo/git_KU/wwt/WWT_PBE/safe')

% % NB this scheme works only for Multivariate Normal Distribution. not other arbitrary multivariate
% distributions


% Numerical solution to the PBE describing cell growth and division. 4 different cases:

% 1) f(x) = g, b(x) = k  ::  GrowthRate='constant', BreakageRate='constant'
% 2) f(x) = g*x, x(b) = k  ::  GrowthRate='linear', BreakageRate='constant'
% 3) f(x) = g, b(x) = k*x  ::  GrowthRate='constant', BreakageRate='linear'
% 4) f(x) = g*x, b(x) = k*x  ::  GrowthRate='linear', BreakageRate='linear'

nD = 1; % number of dimensions of the PBE (= no of distributed particle attributes)
nChar = 1; % number of characteristics 
nCells = 1000; % number of discretiztion points for CAT
GrowthRate = 'linear';%''zero';%'constant'; %'zero'; %'constant';
g = 1; %growth rate, this instead of a vector x1_all
BreakageRate = 'zero';%'zero';%'constant';%'zero';%='constant';
k = 1;%breakage rate


% SCHEME followed
% 1 - MC sampling of non-growth related coordinate x1
% 2 - calculate specific probabilities p and weights w associated to all chosen x1_0 values
% 3 - discretized grid of growth related coordinate x2
% 4 - calculate the joint probability distribution N_0, actually not used
% 5 - calculate the initial marginal distribution N_0~ ti, look to Waldherr hybrid paper
% 6 - calculate trajectories of (x1, w_x1) integrating (dx1dt, dw_x1dt) over time, starting from initial (x1_0, w_x1_0) (MOC, The last characteristic ODE, referring to NDF, is not solved)
% 7 - use (x1, w_x1) to calculate the growth of x2 within the CAT
% 8 - using CAT, come up with N~
% 9 - use N~ to calculate the 2 marginals N_x1 and N_x2, and the moments

% Individual Model
% dxdt = g; linear growth 
% linear breakage rate
% 2 equally sized particles per breakage event
x1Sampling = 'uniform'; % 'MH-sampling', 'normal'
Ntot_0 = 1;

twoPointsFormula = false;
implementGrowth = true;
implementBreakage = true;
thr = 0; % State variables not allowed to go below this threshold
epsylon = 0;
lostFirstMoment = false;
plot3Darmin = true;
[epsyBreak, epsyGrow, nInconsGrowth, nInconsBreak,lostFirstMoment_lb, lostFirstMoment_ub, approx_lb_max, approx_ub_max, tFlaglb, tFlagub] = deal(0);

% growth rate (x1) and cell size (x2) disrtribution parameters
mu_x1 = g;
sigma_x1 = 0.2;
x1Min = 0;
x1Max = 2; 
% --
mu_x2 = 3 ; 
sigma_x2 = 0.5;
x2Min = 0;
x2Max = 800;

if isequal(nD,1)
 nChar = 1;
 x1Min = mu_x1;
 x1Max = mu_x1;
end





[M_Nti_0_full,x2_all,x2_b,x2_lb,x2_ub,x2_x0,N_x0,x2_widths,ghost1,ghostEnd] = defineIC_2D_PBE();
% M_Nti_0_full has each column like [Joint_tilde_i;x1_i;w1_i]. number of columns equals the number of characteristics

% Define breakage kernels for all sizes
[BrRateKer_b,BrRateKer_v] = BrekageKernel(x2_all,x2_lb,x2_ub,nCellsGh);

% time window
tend = 5;
DT = tend/100;
tspan = 0:DT:tend;
tpoints = length(tspan);
step = 0;
M_NisTilde = true; % Used for the ode solution
Nti_res = zeros(tpoints,nCombo);
[x1_all_res, w1_all_res] = deal(zeros(tpoints,nChar));

for aChar = 1:nChar
N_0 = M_Nti_0_full(:,aChar);
% ode system
options = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol', 1e-8); 
[t, N_res] = ode15s(@(t,N_all) odeFun_2D_PBE(t,N_all,x2_all,x2_b,x2_lb,x2_ub,x2_x0,N_x0,x2_widths,BrRateKer_b,BrRateKer_v,ghost1,ghostEnd), tspan, N_0,options);
Nti_res(:,((aChar-1)*(nCellsGh)+(1:nCellsGh))) = N_res(:,(1:nCellsGh));
x1_all_res(:,aChar) = N_res(:,(nCellsGh+1)); 
w1_all_res(:,aChar) = N_res(:,end);
end
Mx2 = zeros(tpoints,nCellsGh); % Marginal for the y variables (growth-related)
Mx1 = zeros(tpoints,nChar); % Marginal for the z variables (non growth-related)
x2_all_res = removeGh(repmat(x2_all',tpoints,1));

[Ntot_res,moments.x1.norm.first, moments.x1.norm.second, moments.x2.norm.first, moments.x2.norm.second] = deal(zeros(tpoints,1));
for jj = 1:tpoints
 jj_s = strcat(['a',num2str(jj)]);
 onetPointM_Nti = reshape(Nti_res(jj,:),nCellsGh,nChar); % for every time point I reconstruct a matrix
 M_Nti_res.(jj_s) = onetPointM_Nti;
 M_N_res.(jj_s) = onetPointM_Nti./w1_all_res(jj,:);
 Marg_x1.(jj_s) = computeMarg(onetPointM_Nti,M_NisTilde,w1_all_res(jj,:),x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
 Marg_x2.(jj_s) = computeMarg(onetPointM_Nti,M_NisTilde,[],[],nChar,'x2');
 Marg_x1_onetPoint_ti = Marg_x1.(jj_s).*w1_all_res(jj,:); % Fundamental, need to return to tilde before the next line
 Ntot_res(jj) = computeMarg(Marg_x1_onetPoint_ti,M_NisTilde,[],[],nChar,'x2');
 
% -- normalize --
 M_N_norm_res.(jj_s) = M_N_res.(jj_s)/Ntot_res(jj); 
 M_Nti_norm_res.(jj_s) = onetPointM_Nti/Ntot_res(jj); 
 M_Nti_norm_res_jj_s_chack2 = M_N_res.(jj_s).*w1_all_res(jj,:); 
 Marg_x1_norm.(jj_s) = computeMarg(M_Nti_norm_res.(jj_s),M_NisTilde,w1_all_res(jj,:),x2_widths,[],'x1'); % with this I obtain the marginal for x1 Marg_x1 (no tilde)
 Marg_x2_norm.(jj_s) = removeGh(computeMarg(M_Nti_norm_res.(jj_s),M_NisTilde,[],[],nChar,'x2'));
 moments.x1.norm.first(jj) = computeMoments(M_Nti_norm_res.(jj_s),M_NisTilde,x2_all,x2_widths,x1_all_res(jj,:),w1_all_res(jj,:),nChar,'x1',1);
 moments.x1.norm.second(jj) = computeMoments(M_Nti_norm_res.(jj_s),M_NisTilde,x2_all,x2_widths,x1_all_res(jj,:),w1_all_res(jj,:),nChar,'x1',2);
 moments.x2.norm.first(jj) = computeMoments(M_Nti_norm_res.(jj_s),M_NisTilde,x2_all,x2_widths,x1_all_res(jj,:),w1_all_res(jj,:),nChar,'x2',1);
 moments.x2.norm.second(jj) = computeMoments(M_Nti_norm_res.(jj_s),M_NisTilde,x2_all,x2_widths,x1_all_res(jj,:),w1_all_res(jj,:),nChar,'x2',2);
end


%%%%%%%%%%%%% print inconsistencies found  %%%%%%%%%%%%%
printInconsistencies()

%%%%%%%%%%%%%%%%%% plots   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% MOMENTS plots and comparison with analytical solution
anal_1m_Res = analyticalMoments(tspan,mu_x2,GrowthRate,BreakageRate);

figure 
plot(tspan,moments.x2.norm.first,'LineWidth', 4)
hold on
plot(tspan,anal_1m_Res,':','LineWidth',4)
plotname = {'First moment \mu_1',strcat(['Growth rate: ',GrowthRate,' Breakage rate: ',BreakageRate])};
% title(plotname,'FontSize', 18);
lgd = legend('Numerical','Analytical');
lgd.FontSize = 18;
xlabel('Time [-]','FontSize', 18) 
ylabel('\mu_1 [-]','FontSize', 25)
cax = gca;
xtick = linspace(cax.XLim(1),cax.XLim(2),5);
ytick = linspace(cax.YLim(1),cax.YLim(2),5);
set(gca, 'XTick',xtick , 'XTickLabel',xtick,'FontSize',18)
set(gca, 'YTick',ytick , 'YTickLabel',ytick,'FontSize',18)
filename = strcat(['M1_GR',GrowthRate,'_BR',BreakageRate]);
saveas(gcf,strcat([filename,'.pdf']))

figure 
plot(tspan,moments.x2.norm.second,'LineWidth', 2)
plotname = {'Second moment',strcat(['Growth rate: ',GrowthRate,' Breakage rate: ',BreakageRate])};
title(plotname)
filename = strcat(['M2_GR',GrowthRate,'_BR',BreakageRate]);
saveas(gcf,strcat([filename,'.pdf']))

figure
plot(tspan,Ntot_res,'LineWidth', 2) 
plotname = strcat(['Ntot ',' [GrowthRate:',GrowthRate,' BrRate:',BreakageRate,']']);
title(plotname)
filename = strcat(['Ntot_GR',GrowthRate,'_BR',BreakageRate]);
saveas(gcf,strcat([filename,'.pdf']))
return

function [] = plotMovie(tpoints,x1_all_res,x2_all_res,Marg_x1_norm,Marg_x2_norm)
figure
for jj = 1:tpoints
 jj_s = strcat(['a',num2str(jj)]);
 plot(x1_all_res(jj,:), Marg_x1_norm.(jj_s), 'LineWidth', 2) 
% hold on
 plot(x2_all_res(jj,:), Marg_x2_norm.(jj_s), 'LineWidth', 2)  
% hold off
 legend('Marg x1')
 axis square
 title(jj_s)
 drawnow;
 pause(0.25);
end
return



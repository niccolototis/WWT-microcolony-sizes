function [M_N_TI_0_full,P] = defineIC(aSt,P)
global speciesPBE adhDistribution nChar  initialNDF   nCells  spNames   
global impulseDistribution begN endN  importanceSampling speciesUnif distrib_at_startupScenario  nSp myPlotOff  size_N_t0_aSpPBE_  
% [1] - Steffen waldherr hybryd simulations in heterogeneous cell
% populations
% [2] - Physically Based Rendering: From Theory to Implementation
% ---------------------------------------
% Please note, in the following functions: 
% 1) x1 always refer to the concentration of adhesins ADH
% 2) ! NB mu_v can be either a volume or biomass, depending of ModelBiomass = true; in Main file, and the
% assigment in defineParams_WWT_HETUnif.m
% ---------------------------------------
% NB computeMarg is used to compute the marginal for x1 only! Instead
% integration over x2 coordinate is done just by summming up, because the
% vector N_data already contains the zeroth moments for all
% subdomains
% ---------------------------------------
% NB here the x2-coordinates have already been RESCALED so that they are
% between 0 and 1 = > x2_bar_RESC0to1 has no real physical sense, I should multiply
% it by the range of original coordinates to obtain the real x2bar



size_N_t0_aSpPBE_ = (nCells+2); % +2 is for the w_i and x_i
[begN,endN] = deal(zeros(nSp,1));

for aSpPBE = speciesPBE
 sp = char(spNames{aSpPBE});
 begN(aSpPBE) = speciesUnif + (find(speciesPBE==aSpPBE)-1)*size_N_t0_aSpPBE_ +1;     
 endN(aSpPBE) = begN(aSpPBE)+size_N_t0_aSpPBE_-1;
 
 [mu_x1,x1Min,x1Max,mu_x2,x2_b,x2_lb,x2_ub,x2_all,x2_widths,mu,sigma,N_data] = assignStruc(aSt.(sp));

 if aSt.(sp).mu_x2<(x2_all(1)+x2_lb(1))/2 % they should sum up to one aprt from where I have no births
    ME = MException('pdf will be all zeros: NEED smaller cells close to the origin >> increase diamInitPart in defineParams_WWT.m, reduce diamMax or increase the number of cells nCells');
    throw(ME)
 end


 %% 2 - use Metrppolis Hastings algorithm for MC sampling from the marginal distribution (of the joint) for non-growth related coordinate x1 (integrating over all x2_all points)
 % smpl = mhsample(start,nsamples,'pdf',pdf,'proppdf',proppdf, 'proprnd',proprnd)
 if importanceSampling
        start = mu(1);  %  % IMPO DO NOT USE mu_adh because in the case of lognormal distribution would not work!
        discarded = 200;
        nsamples = nChar;
        % f{1} = @(x) x;
        nx1 = 1;
        MargForWhat = 'x1';
        % M_N_0 = computeJoint(x1pts,x2pts,nx1,nx2, mu, sigma);
        % [Marg] = computeMarg(M_N_0,w1_all,x2_widths,nx1,MargForWhat)
        thisPdf = @(x1pt) computeMarg(computeJoint(x1pt,x2_all,nx1,nCells,mu, sigma),[],x2_widths,nx1,MargForWhat); %
        delta = sigma(1); % IMPO DO NOT USE sigma_adh because in the case of lognormal distribution it causes error!
        proppdf = @(x,y) normpdf(x,y,delta); %The proposal distribution q(x,y) gives the probability density for choosing x as the next point when y is the current point. Use here current point y as the mean of the proposal distrib function == normpdf
        proprnd = @(y) y + randn*delta;  % proprnd defines the random number generator for the proposal distribution (qui genero x a seconda di y
        [x1_all_sp_norm,~] = mhsample(start,nsamples,'pdf',thisPdf,'proppdf',proppdf, 'proprnd',proprnd,'symmetric',false,'burnin',discarded);
        x1_all_sp_norm = sort(x1_all_sp_norm);
        switch adhDistribution
            case 'normal'
                x1_all = x1_all_sp_norm;
            case 'lognormal'
                x1_all = exp(x1_all_sp_norm);
        end
        
 if ~myPlotOff
    figure;
    h = histfit(x1_all,50,'kernel'); % histfit(b,10,'kernel')
    hold on
    switch adhDistribution
        case 'normal'
             plot(x1_all,nsamples/20*normpdf(x1_all,aSt.(sp).mu_x1,aSt.(sp).sigma_x1),'LineWidth',2);
        case 'lognormal'
%              plot(x1_all,nsamples/20*lognpdf(x1_all,mu_adh_NfromLog.(sp),sigma_adh_NfromLog.(sp)),'LineWidth',2);
    end
%             plot(x1_all,pdf(x1_all),'g','LineWidth',2);
    h(1).FaceColor = [.8 .8 1];
    hold off
 end
 else
    x1_all = (x1Min:(x1Max/(nChar-1)):x1Max);
 end
 aSt.(sp).x1_all = x1_all;
 
 %% 3 - Calculate the initial joint probability density function
 % To create the initial ndf N I first define MVNPDF Multivariate normal
 % probability density function (pdf), and then I multiply it by the
 % total number of particles
 % NB: the normalization (optional) is lated done on the Ñ, not on
 % N
 switch initialNDF
    case 'normal_2D'
        P.M_N_0 = computeJoint(aSt.(sp));  %IMPO! this x1_all_sp_norm needs to be here
    case 'impulse_2D'
        if impulseDistribution % replace M_N_0 witht the impulse function
            P.M_N_0 = x2_all*x1_all*0;
            [minVal,irow] = min(abs(x2_all-mu_x2));
            [minVal,icol] = min(abs(x1_all-mu_x1));
            P.M_N_0(irow,:) = (1/P.nChar)/x2_widths(irow); 
        end
    case 'x2MarginalByKDEstimate'
        if strcmp(P.IC_at,'startupScenario')
            switch distrib_at_startupScenario
                case 'impulse'
                    P.M_N_0 = zeros(size(x2_all));
                    P.M_N_0(1) = 1; % I am putting 
                case 'lognormal'
                    P.M_N_0 = computeJoint(aSt.(sp),'lognorm');  %IMPO! this x1_all_sp_norm needs to be here
                case 'normal'
                    P.M_N_0 = computeJoint(aSt.(sp),'norm');  %IMPO! this x1_all_sp_norm needs to be here
            end
        else
            P.M_N_0 = N_data;
            reallyEqual(sum(P.M_N_0),1,'should be 1')
        end
 end 
      
    
 if length(P.M_N_0(P.M_N_0~= 0))==0 % if it is all zeros
    ME = MException(strjon(['The initial joint distribution has not been defined' ... 
    'correctly, review computeJoint()']));
    throw(ME)
 end
 
 %% 4 - Calculate the initial marginal frequency for the chosen x1_all, calculate associated weights 
 MargForWhat = 'x1'; % function [Marg]=computeMarg(M_N_0Mx1_0.(sp),w1_all,x2_widths,nx1,MargForWhat)
 Mx1_0 = sum(P.M_N_0,1);
 
 %
 int_Joint_NDFnorm_0 = 1;
 w1_all = (int_Joint_NDFnorm_0./Mx1_0); % point 1) right column of the second page of [1]   
 
 %{
 Following [2]
 I = integral to be estimated with Monte Carlo method
 I = 1/N*symsum( f(x(i))/p(x(i)) ),i,0,N)
 p(x(i)) = c*f(x(i)) % I want the probability of my samples to have more or less the same distribution of
 the pdf I want to calculate the integral of
 by normalization c = 1/int(f(x)0
 Since c is a constant, each estimate has the same value, and the variance is zero! Of course, this is ludicrous since we wouldn?t bother using Monte Carlo if we could integrate f(x) directly. However, if a density p(x) can be found that is similar in shape to f(x), variance decreases.
 %}

 
 
 %% 5 - Define ndf(t = 0)
 % I consider that at time t = 0 joint_ndf==joint_pdf, with zeroth moment=1. From the first moment \bar{x2}=E[x2] 
 % I compute the initial number of particles joyFigs, which is used
 % a posteriori as a rescaling factor to come up with the real ndf(t)
 
 % - calculate Ñ as in equation 10 in [1]
 P.M_N_TI_0 = P.M_N_0.*w1_all; % Each column multiplied by a specific weight

 % By computing the first moment I get the average particle size x2_bar_RESC0to1, 
 % then I get Ntot by InitSize_all/x2_bar_RESC0to1. joyFigs is later used in 
 % the main after the PBE solution is computed 
 
 P.M_NisTilde = true;  
 % - Integral along the rows (x2 coord)
 firstMom_x2 = sum(P.M_N_TI_0.*x2_all,1);
 % - Integral across the columns (x1 coord)
 x2_bar_RESC0to1 = computeMarg(firstMom_x2,P.M_NisTilde,[],[],P.nChar,'x2');
 P.N_tot_0_all.(sp) = P.InitSize_all(aSpPBE)/x2_bar_RESC0to1;
 
 
 
 %% 7 - Check that the integrals are implemented correctly and the zeroth and first moments matches Ntot and InitSize_all 
 int_Joint_NDF_0 = int_Joint_NDFnorm_0* P.N_tot_0_all.(sp); 
 checkIntegralCorrectness(P.M_N_TI_0,P.nChar,int_Joint_NDFnorm_0);
 M_N_TI_0_ndf = P.M_N_TI_0*P.N_tot_0_all.(sp);
 checkIntegralCorrectness(M_N_TI_0_ndf,P.nChar,int_Joint_NDF_0);
 

 
 %% 7 - Create output    
 if P.normalized_N_TI_t0
    M_N_TI_0_full.(sp) = [P.M_N_TI_0; x1_all; w1_all];  
    disp('At t = 0 the initial ndf is defined as a probability density function')
 else
    M_N_TI_0_full.(sp) = [M_N_TI_0_ndf; x1_all; w1_all];  
    disp('At t = 0 the initial ndf is NOT a probability density function')
 end
 plotIC(M_N_TI_0_ndf,aSt.(sp),x2_bar_RESC0to1)
end
return 

function [mu_x1,x1Min,x1Max,mu_x2,x2_b,x2_lb,x2_ub,x2_all,x2_widths,mu,sigma,N_data] = assignStruc(aStp)
global initialNDF IC_at
 mu_x1 = aStp.mu_x1;
 x1Min = aStp.x1Min;
 x1Max = aStp.x1Max;
 mu_x2 = aStp.mu_x2;
 x2_b = aStp.x2_b;
 x2_lb = aStp.x2_lb;
 x2_ub = aStp.x2_ub;
 x2_all = aStp.x2_all;
 x2_widths = aStp.x2_widths;
 mu = aStp.mu;
 sigma = aStp.sigma;
 if strcmp(initialNDF,'x2MarginalByKDEstimate') && ~strcmp(IC_at,'startupScenario')
    N_data = aStp.N_data;
 else
    N_data = [];
 end
return

function [] = checkIntegralCorrectness(M_N_TI_0,nChar,zeroth_TRUE)
 % This check is used both when M_N_0 given as input is the ndf or the
 % normalized_ndf(which is a pdf)
 % compare the zeroth moments computed in different ways
 Marg_x1_TI = sum(M_N_TI_0,1);
 M_NisTilde = true;
% Marg_x1_TI = Marg_x1.*w1_all; % Fundamental, need to return to tilde before the next line
 zeroth_1 = computeMarg(Marg_x1_TI,M_NisTilde,[],[],nChar,'x2');
 % 
 Marg_x2 = computeMarg(M_N_TI_0,M_NisTilde,[],[],nChar,'x2');
 zeroth_2 = sum(Marg_x2,1);
 % 
 
 reallyEqual(zeroth_1,zeroth_TRUE,'These zeroth moments ');
 reallyEqual(zeroth_2,zeroth_TRUE,'These zeroth moments ');
return

function [] = plotIC(M_N_TI_0,aStp,x2_bar_RESC0to1)
global speciesPBE spNames scaleDots   myPlotOff  nChar
% Marginal x2
MargForWhat = 'x2';
% MC integration of marginal for x_1 Mx1_0 http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Importance_Sampling.html
M_NisTilde = true;
Mx2_0 = computeMarg(M_N_TI_0,M_NisTilde,[],aStp.x2_widths,nChar,MargForWhat);

if ~myPlotOff
figure
for aSpPBE = speciesPBE
 sp = char(spNames{aSpPBE});
 plot(aStp.x2_all, Mx2_0, '-', 'LineWidth', 2)
 hold on      
end
xline(x2_bar_RESC0to1,'r-');
xline(aStp.x2SCell,'b-');
plot(aStp.x2_all, 0, scaleDots, 'LineWidth', 2)
% plot(aStp.x2_lb(2),0,'rx'); % plot range
% plot(aStp.x2_ub(end-1),0,'rx'); 
hold off
title('Marginal for v')
legend(spNames(speciesPBE));
end

% Marginal x1
Mx1_0 = sum(M_N_TI_0,1);
 
if ~myPlotOff 
figure
for aSpPBE = speciesPBE
 sp = char(spNames{aSpPBE});
 plot(aStp.x1_all, Mx1_0, '-', 'LineWidth', 2)
 hold on      
end
plot(aStp.x1_all, 0, scaleDots, 'LineWidth', 2) 
hold off
title('Marginal for adh')
legend(spNames(speciesPBE));
end
return


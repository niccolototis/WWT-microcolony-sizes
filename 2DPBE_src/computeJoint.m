function [M_N_0] = computeJoint(aStp,whatdistr) %x1pts,x2pts,nx1,nx2, mu, sigma)
global nChar 
 nx1 = nChar;
 x1pts = aStp.x1_all;
 x2_b = aStp.x2_b;
 x2_lb = aStp.x2_lb;
 x2_ub = aStp.x2_ub;
 nx2_b = length(x2_b);
 
 % - with meshgrid in X I get one x1pts per column, in Y I get one x2_b per row   
 [X, Y] = meshgrid(x1pts,x2_b);
 jointCoord = reshape([X, Y],nx1*nx2_b,2);
 if nChar ==1
    % I multiply all coordinates and statistics by this factor
    % otherwise with cdf or pdf I get only zeros (apparently because
    % the with of each interval is too small?)
    help = 1000;
    switch whatdistr
        case 'norm'
            N_x2_b = normcdf(jointCoord(:,2)*help, aStp.mu_x2*help, aStp.sigma_x2*help);
        case 'lognorm'
            N_x2_b = logncdf(jointCoord(:,2)*help, aStp.mu_x2*help, aStp.sigma_x2*help);
    end
 else
    error('Only the case of nChar ==1 has been implemented for this feature');
 end
 N_at_x2_all_toBeRescaled = N_x2_b(2:end)-N_x2_b(1:end-1);
 % I rescale it because I need it to sum up to 1
 M_N_0 = N_at_x2_all_toBeRescaled/sum(N_at_x2_all_toBeRescaled); 
return

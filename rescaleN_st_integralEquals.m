function [N_at_x2_all,proportionalityK]=rescaleN_st_integralEquals(N_at_x2_all,mustintegrateTo)
    integralOfTheMarg=sum(N_at_x2_all,1);
    proportionalityK=(mustintegrateTo./integralOfTheMarg)';
    % - Both N_at_x2_all and proportionalityK need to be in column
    N_at_x2_all=N_at_x2_all*proportionalityK;
return

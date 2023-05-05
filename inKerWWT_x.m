function outp = inKerWWT_x(x,BorV)
global redistrFn 
% Considering daughter particle (volume vD) produced by the breakage of a mother particle (volume vM), with vD<= vM 
% intBrKerWWT_V2 calculates:
% 1) the zeroth moment (B) = the probability of obtaining a particle of vD
% in [vD-,vD+] from a particle of vM
% 2) the first moment, i.e. the expected value of vD for a particle generated in [vD-,vD+]
% To model erosion, the arcsine distribution is used:
% pdf = 1./(pi*sqrt(x.*(1-x)))
% cdf = 2/pi*asin(sqrt(x))

switch redistrFn
 case 'U-shape' 
    switch BorV           
    case 'B' 
        % look to compareIntegrals.m for reference:
        m0_x = (2/pi).*asin(sqrt(x)); 
        outp = m0_x;
    case 'V' 
        % look to compareIntegrals.m for reference
        m1_x = 1/2-(sqrt(x-x.^2)+asin(sqrt(1-x)))./pi;
        outp = m1_x;
    end
end
return 



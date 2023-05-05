function outp = intBrKerWWT(v_x,BorV,vDmax,vMin)%(vD,vDmax,aSp) %vD=vDmax=parameter %breakage kernel [or daughter size distribution] ?(probability density function for the formation of particles of size vD from particles of size vDmax). Here time independent
% thi is how I am calling it: intBrKerWWT(ubs,'B',vDmax,lbs(1))
global redistrFn FirstRefNoBreak
% intBrker calculates the integral of the probability of obtaining a daughter particle (with volume vD)
% from the breakage of a mother particle (with volume vDmax). The integral is calculated for each ALLOWED
% vD cell sizes (vD<vDmax by assumption)
%%% Breakage kernel (resistribution function)
% Arcsine distribution
% pdf = 1./(pi*sqrt(x.*(1-x)));
% cdf = 2/pi*asin(sqrt(x)); % this is the analythical solution of the integral. I want to evaluate this at
% each boundary
 vDinterval = vDmax-vMin;
 if vDinterval ==0
    if ~FirstRefNoBreak
        ME = MException('This should happen only when the size of the mother particle coincides with the lowest possible daugter size');
        throw(ME)            
    end 
    outp = 0;
 else
 switch redistrFn
    case 'constant'
        switch BorV             
        case 'B' 
            outp = v_x./vDinterval;
        case 'V' 
            outp = (v_x.^2)/(2*vDinterval); %cumulative volume = birts*expected volume. Expected volume is arithmetic mean [(vB+vA)/2] , that I have to multiply by births *(vB-vA)/vDmax . .. secondo termine e' la probabilita. if the probability is constant over the spaces daughter cell volume vD, then in the volume increase due by birth in cell_i has to equal the arithmetic mean of volume in cell_1 * Bb_i
        end
    case 'U-shape' 
        v_xN = (v_x-vMin)./vDinterval;% normalize daughter size range between 0 and 1 (outside this range the integral returns imaginary numbers)
        if ~isempty(v_xN(v_xN>1|v_xN<0))>0 
            ME = MException('should be normalized this kernel fcn is defined only between 0 and 1!');
            throw(ME)            
        end         
        switch BorV           
        case 'B' 
            outp = 2/pi*asin(sqrt(v_xN));  % Output is a probability normalization does not change anything here 
        case 'V' 
            % @@NT impo here I have to use v_v_xN normalized BUT also v_x here
            exp_v_xN = 1/2-(sqrt(v_xN - v_xN.^2) + acos(sqrt(v_xN)))./pi; % expected value of the normalized v_xN 
            exp_v_x = exp_v_xN*vDinterval+vMin;
            outp = exp_v_x; % go back from the normalized expected value to the original expected value. should work because integration is a linear operator!                  
        end
 end
 end
return 

function retVal = grRate(x2_all,g)
global GrowthRate
 switch GrowthRate
    case 'zero'
        L = 0;
        g = 0;
    case 'constant'
        L = 0;     
    case 'linear'
        L = 1;
    case 'quadratic'
        L = 2;
    case 'cubic'
        L = 3;
 end
 retVal = g * x2_all.^L;  
return



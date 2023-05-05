function retVal = BrRate(v,aSp)
global BreakageRate k
 switch BreakageRate
    case 'zero'
        L = 0;
        k = 0;
    case 'constant'
        L = 0;     
    case 'linear'
        L = 1;
    case 'quadratic'
        L = 2;
    case 'cubic'
        L = 3;
 end
 retVal = k * v.^L;  
return

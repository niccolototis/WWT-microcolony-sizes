function retVal = BrRateWWT(v,P)
 switch P.BreakageRate
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
 retVal = P.lambda/P.adhesins .* v.^L;  
return


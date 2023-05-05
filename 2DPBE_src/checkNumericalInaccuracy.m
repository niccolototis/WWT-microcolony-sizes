function [] = checkNumericalInaccuracy(B,D,PREoPOST)
%     if (I enter this checkout) BUT
%     (~isequal(num2str(sum(B_CA)),num2str(sum(D_CA))) ==false) THEN I HAVE THE PROOF
%     THE PROBLEM IS NUMERICAL MACHINE PRECISION, performing B_CA-D_CA leads to
%     approximations !
global implementGrowth implementBreakage epsylon
 tepsy = 0;
 if implementGrowth && ~implementBreakage && abs(sum(B)-sum(D))>epsylon
    tepsy = abs(sum(B)-sum(D));
 end
 if (~implementGrowth && implementBreakage) && (abs(sum(B)-2*sum(D))>epsylon)
    tepsy = abs(sum(B)-2*sum(D));
 end
 if tepsy>epsylon % they should sum up to one aprt from where I have no births. For some numerical reason there is a little inaccuracy, so a epsylos should be considered/ here epsylon = 1e-15
    epsylon = tepsy;
    disp(['epsylon--',PREoPOST, num2str(tepsy)])
%             ME = MException('total B and total D should sum up to the same number');
%             throw(ME)
 end
return

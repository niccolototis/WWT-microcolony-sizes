function N = fixLowerThreshold(N,thr)
global minNeg
 if ~isempty(N(N<thr))
    minNegNew = min(N(N<thr));
    if minNegNew<minNeg
        disp(['min NEGATIVE ', num2str(minNegNew)])
        minNeg = minNegNew;
    end
    N(N<thr) = thr;
 end
return

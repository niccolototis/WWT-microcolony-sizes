function N = fixLowerThreshold_fast(N,P)
global minNeg step
if ~isempty(N(N<P.thr))
 minNegNew = min(N(N<P.thr));
 indMin = find(N==minNegNew);
 if strcmp(P.isRunning,'PBE')
    if ~P.identify_lambda && ~P.sweep_lambda
        if minNegNew<minNeg
            disp(['min NEGATIVE ', num2str(minNegNew),' position ',num2str(find(N ==minNegNew)), ' step=',num2str(step)])
            minNeg = minNegNew;
        end
    end
 end    
 tooMuch = -0.1;
 if minNegNew<tooMuch
    ME = MException(strcat(['Pay attention, values going negative below ',num2str(tooMuch),'. Either 1) there '...
    'is something going wrong; 2) nothing is going wrong, it is just a glitch of the odesolver that is making some punctual '...
    'attempts of integration with very different values. When opt 1 can be excluded, you can comment out these lines']));
    throw(ME);
 end
 N(N<P.thr) = P.thr;
end
return

function P = paramsToModify(P)
global parametersEffectOnMSD noTrials 
 if parametersEffectOnMSD
    P.impro_chosen.QIN=[210,210,210]; 
    P.impro_chosen.SRT=[15,10,5]; % Only SRT is set to change
    P.impro_chosen.O2_c_ref=[1.5,1.5,1.5]; 
    P.impro_chosen.cycleVariant=[0 0 0];
    fN = fieldnames(P.impro_chosen);
    for q = 1:length(fN)
        if q>1
            if ~isequal(ll,length(P.impro_chosen.(fN{q})))
                error('should have all the same lenght')
            else
                ll = length(P.impro_chosen.(fN{q}));
            end
        else
            ll = length(P.impro_chosen.(fN{q}));
        end
    end

    % update noTrials if I am choosing it manually it needs to be updated
    noTrials = ll;
 end
end

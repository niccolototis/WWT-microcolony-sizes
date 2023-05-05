function results = collectResultOnePartitionPBE(first,last,aSt,P,M_Nti_0_full,t_all_ASM1seq,S_c_V_whole_tall,Ker_B,Ker_V,sp,Q)
% 6 different values of the objective function, computed with the different metrics and criteria
 results = zeros(last-first,6); 
 for ii = first:last
    if P.debugParallel
        P.lambda = P.lambdasToSweep(ii);
        aSt = runPBEforOneParam(aSt,P,M_Nti_0_full,t_all_ASM1seq,S_c_V_whole_tall,Ker_B,Ker_V,sp);
        oneResult = [];
    else
        try
            P.lambda = P.lambdasToSweep(ii);
            aSt = runPBEforOneParam(aSt,P,M_Nti_0_full,t_all_ASM1seq,S_c_V_whole_tall,Ker_B,Ker_V,sp);
            oneResult = [];
        catch
            warning('Problem during integration with parameters.  Assigning a value of NaN to the output of the run.');
            oneResult = NaN;
        end
        % Plotting the final value of HET biomass on the z axis of thisFig
        send(Q,[ii,oneResult]);
    end
    % Depending on what I want to return I initialize the matrix/vector
    % with output results to be fetched
    results(ii-first+1,:) = oneResult;
 end
return

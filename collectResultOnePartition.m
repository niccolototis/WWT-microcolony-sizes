function results = collectResultOnePartition(first,last,P,Q)
 % Define Helper Functions
 % Define a helper function that solves the Lorenz system on a partition of the parameters to explore. Send intermediate results to the MATLAB client by using the send function on the DataQueue object.
 results = zeros(last-first,1);
 switch P.figure
    case 'surface'
        results = zeros(last-first,1);
    case 'flat'
        tPts = P.ntpts_ASM1seq;
        if strcmp(P.plot,'both_HET_AUT')
            tPts = tPts*2;
        end
        results = zeros(last-first,tPts);
 end
 for ii = first:last-1
    try
        oneResult = runASM1seqForOneParam(ii,P)';
    catch
        warning('Problem during integration with parameters.  Assigning a value of NaN to the output of the run.');
        oneResult = NaN;
    end
    %plotting the final value of HET biomass on the z axis of the
    %thisFig
    send(Q,[ii,oneResult]);
    results(ii-first+1,:) = oneResult;
 end
return

function [] = computeBestLambda_fn(P,for_or_parfor)
global IC_at folderOutput cycleModified parametersEffectOnMSD cycleModStr improvStr lambdasToSweep howDistancesComputed nDecimalsToIncludeNDFsampling
switch IC_at
 case 'SS_computed'
    lambdasFiles = dir([folderOutput,'P_allLambdas/']);
    lambdasFiles = {lambdasFiles.name};
    lambdasFiles = lambdasFiles(contains(lambdasFiles,'P_lambda'));
    if cycleModified || parametersEffectOnMSD
        if cycleModified
            lambdasFiles = lambdasFiles(contains(lambdasFiles,cycleModStr));
        end
        if parametersEffectOnMSD
            lambdasFiles = lambdasFiles(contains(lambdasFiles,improvStr));
        end
    else
        lambdasFiles = lambdasFiles(~contains(lambdasFiles,'_cycle_'));
    end
    if ~isempty(lambdasFiles)
        lambdasFiles = lambdasFiles{1};
    else
        error('Results are missing')                
    end

    % This just to load back P.lambdasToSweep
    load(strcat([folderOutput,'P_allLambdas/' lambdasFiles]),'P');
    P = PBEsettings(P);
    nCells_ = P.nCells; T_sim_=P.T_sim;
    % Criterio con cui l'ottimo viene scelto, considerando quale param
    for i = 1:length(lambdasToSweep)    
        Lambda_i = lambdasToSweep(i);
        labLambda = ['lambda_' num2str(Lambda_i)]; 
        P = [];
        load(strcat([folderOutput,'P_allLambdas/P_',labLambda,improvStr,cycleModStr,'.mat']),'P');
        filename_aSt = strcat([folderOutput,'aSt_allLambdas/aSt_',labLambda,'_T_sim',num2str(T_sim_),'_nCells',num2str(nCells_),improvStr,cycleModStr,'.mat']);
        load(filename_aSt,'aSt');
        % Compute the 3 statistics        
        labLambda = strrep(labLambda,'.','_');
        switch howDistancesComputed
            case 'tptsApart'
                Dist = computeSimilarityTest_tptsApart(aSt,Lambda_i);
            case 'tptsAllin'
                switch for_or_parfor
                    case 'for'
                        Dist = computeSimilarityTest_tptsAllin_for(aSt,Lambda_i);
                    case 'parfor'
                        Dist = computeSimilarityTest_tptsAllin_parfor(aSt,Lambda_i);
                end
        end
        Dist_lambdasToSweep.(labLambda) = Dist;
        save([folderOutput 'Dist_lambdasToSweep_' howDistancesComputed '_' num2str(nDecimalsToIncludeNDFsampling) 'dec' improvStr cycleModStr '.mat'],'Dist_lambdasToSweep')
    end    
 case 'startupScenario'
    labLambda = ['Lambda_' num2str(optimal_lambda)];
    load(strcat([P.folderOutput,'P_allLambdas/P_', distrib_at_startupScenario ,'_',labLambda,improvStr,cycleModStr,'.mat']),'P'); 
    P = PBEsettings(P);
    % - load PBE data
    aSt = [];
    load(strcat([P.folderOutput,'aSt_allLambdas/aSt_', distrib_at_startupScenario ,'_',labLambda,'_T_sim',num2str(P.T_sim),'_nCells',num2str(P.nCells),improvStr,cycleModStr,'.mat']),'aSt');
    % Compute the 3 statistics        
    labLambda = strrep(labLambda,'.','_');
    switch howDistancesComputed
        case 'tptsApart'
            Dist = computeSimilarityTest_tptsApart(aSt,optimal_lambda);
        case 'tptsAllin'
            switch for_or_parfor
                case 'for'
                    Dist = computeSimilarityTest_tptsAllin_for(aSt,Lambda_i);
                case 'parfor'
                    Dist = computeSimilarityTest_tptsAllin_parfor(aSt,Lambda_i);
            end
    end
    Dist_lambdasToSweep.(labLambda) = Dist;
    save([folderOutput 'Dist_lambdasToSweep_' howDistancesComputed '_' num2str(nDecimalsToIncludeNDFsampling) 'dec' improvStr cycleModStr '.mat'],'Dist_lambdasToSweep')
end
end

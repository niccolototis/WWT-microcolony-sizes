function [P,varargout] = parallelRuns_ASM1seq(P)
% Create a figure object, and set 'Visible' to true so that it opens in a new window, outside of the live script. To visualize the results of the parameter sweep, create a thisFig plot. Note that initializing the Z component of the thisFig with NaN creates an empty plot.
 P.runningParallel = true;
 P.nRuns_tot = P.nRuns;
 while P.nRuns_tot>0       
    P.nRuns = min(P.nRuns_tot,P.batchSizeSweep);
    P.nRuns_tot = P.nRuns_tot-P.nRuns;
    switch P.parSetup
        case 'paramSweep'
            switch P.figure
                case 'surface' 
                    thisFig = surf(P.p1,P.p2,NaN(size(P.p1)));
                    xlabel('iXB')
                    ylabel('\mu A','Interpreter','Tex')
                    zlabel('net AUT growth %')
                    P.nRuns = numel(P.p1);   
                case 'flat'
                    P.lgdEntries = [];
                    P.lgdfullNames = [];
            end
        case 'EEMorris'
            P.figMorris = surf((1:(P.EE.k+1)),(1:P.EE.nTraj),NaN(P.EE.nTraj,(P.EE.k+1))); 
            % note that the matrix for the Z coord needs to be given with dimensions in the opposite order
            % e.g: surf(1:5,1:2,rand(2,5))
            xlabel('paramLevel')
            ylabel('trajectory')
            zlabel('relative AUT growth')
            P.title = 'Morris (k+1) runs for each trajectory';
            title(P.title)
            P.nRuns = numel(P.EE.X_Params(:,1,:)); 
            % - Prepare matrix of used parameters
            ii_all = 1:P.nRuns;
            % The P.EE.k+1 x P.EE.nTraj face of the P.EE.X_Params
            % cube isc scrolled "vertically" going down the columns
            % (column-wise)
            %: if I do this than I can use the normal reshape (column-wisw)
            [ind_x,ind_z] = ind2sub([P.EE.k+1 P.EE.nTraj],ii_all);
            % - check that the order is consistent with the reshape
            % that I do later with fetchoutputs, so that the
            % correct relationships between (k,nTraj) indexes and
            % the index inside fetchedOutputs is maintained
            % check_order = string(ind_x(:))+'_'+string(ind_z(:));
            % disp(reshape(check_order,P.EE.k+1,P.EE.nTraj))  
            usedParams = zeros(P.EE.k,P.nRuns); % in rows I have all the params, in colums the different sets used
            for ii = ii_all
                usedParams(:,ii) = P.EE.X_Params(ind_x(ii),:,ind_z(ii));
            end
            nn = usedParams';
            [u,I,~] = unique(nn, 'rows', 'first');
            hasDuplicates = size(u,1) < size(nn,1);
            ixDupRows = setdiff(1:size(nn,1), I);
            dupRowValues = nn(ixDupRows,:);
            if ~isempty(dupRowValues) || ~isempty(ixDupRows) || hasDuplicates
                ME = MException('All columns of usedParams should be different');
                throw(ME);
            end             
    end

    switch P.parallelMode
        case 'parfeval'
            % Perform Parallel Parameter Sweep
        %     refer to https://nl.mathworks.com/help/parallel-computing/run-code-on-parallel-pools.html#:~:text = What%20Is%20a%20Parallel%20Pool,cluster%20in%20your%20parallel%20preferences.

            % After you define the parameters, you can perform the parallel parameter sweep.
            % parfeval works more efficiently when you distribute the workload. To distribute the workload, group the parameters to explore into partitions. For this example, split into uniform partitions of size step by using the colon operator (:). The resulting array partitions contains the boundaries of the partitions. Note that you must add the end point of the last partition.

            tryThis = 1; % or 2
            stepPart = max(1,round(10^(numel(num2str(P.nRuns))-tryThis)));
            partitions = [1:stepPart:P.nRuns, P.nRuns+1]; % numel is the number of array elements

            % For best performance, try to split into partitions that are: 
            % Large enough that the computation time is large compared to the overhead of scheduling the partition.
            % Small enough that there are enough partitions to keep all workers busy.
            % To represent function executions on parallel workers and hold their results, use future objects.

            % - numout the number of output remains 1 even if I am returning a
            % vector from runASM1seqforOneParam
            numout = 1;
            f(1:numel(partitions)-1,numout) = parallel.FevalFuture;

            % Offload computations to parallel workers by using the parfeval function. parallelRuns_ASM1seq is a helper function defined at the end of this script that solves the Lorenz system on a partition of the parameters to explore. It has one output argument, so you must specify 1 as the number of outputs in parfeval.

                % Set Up Parallel Environment
            % Create a pool of parallel workers by using the parpool function.

            if ~isempty(gcp('nocreate'))
                delete(gcp('nocreate'))
            end
            parpool([1 36])

            % To send data from the workers, create a DataQueue object. Set up a function that updates the thisFig plot each time a worker sends data by using the afterEach function. The updateASM1seqPlot function is a supporting function defined at the end of the example.

            Q = parallel.pool.DataQueue;

            afterEach(Q,@(data) updateASM1seqPlot(data,P));
        %     minFuture = afterAll(entry, @(datalgd) mergedLgd(datalgd), 1);

            for ii = 1:numel(partitions)-1
                first = partitions(ii);
                last = partitions(ii+1);
                    f(ii,:) = parfeval(@collectResultOnePartition,1,first,last,P,Q);
            end

            % parfeval does not block MATLAB, so you can continue working while computations take place. The workers compute in parallel and send intermediate results through the DataQueue as soon as they become available.
            % If you want to block MATLAB until parfeval completes, use the wait function on the future objects. Using the wait function is useful when subsequent code depends on the completion of parfeval.
            wait(f);
        case 'parfor'
            tPts = P.ntpts_ASM1seq;
            if strcmp(P.plot,'both_HET_AUT')
                tPts = tPts*2;
            end
            results = zeros(P.nRuns,tPts);
            grNet_AUT_allRuns = zeros(P.nRuns,P.ntpts_ASM1seq);
            X_c_end_allRuns = zeros(P.nRuns,P.nPart);
            S_c_end_allRuns = zeros(P.nRuns,P.nSolub);
            parfor aRun = 1:P.nRuns
                if P.tryCatch                    
                    try
                        [oneResult,~,S_c_tall_gc,X_c_end,S_c_end] = runASM1seqForOneParam(aRun,P);
                    catch
                        warning('Problem during integration with parameters.  Assigning a value of NaN to the output of the run.');
                        [oneResult,X_c_end,S_c_end] = deal(NaN);
                        S_c_tall_gc = nan(1,P.nSolub);
                    end
                else
                    [oneResult,~,S_c_tall_gc,X_c_end,S_c_end] = runASM1seqForOneParam(aRun,P);
                end
                updateASM1seqPlot([aRun,oneResult],P);
                                            
                results(aRun,:) = oneResult;
                grNet_AUT_allRuns(aRun,:) = P.MuA.*(S_c_tall_gc(:,P.sNH4)./(P.KNH+S_c_tall_gc(:,P.sNH4))).*(S_c_tall_gc(:,P.sO2)./(P.KOH+S_c_tall_gc(:,P.sO2)))-P.bA;
                X_c_end_allRuns(aRun,:) = X_c_end;
                S_c_end_allRuns(aRun,:) = S_c_end;
            end  
    end
    delete(gcp('nocreate'))
    P.runningParallel = false;

    switch P.parSetup
        case 'paramSweep'
            varargout{1} = results;
            varargout{2} = grNet_AUT_allRuns;
            varargout{3} = X_c_end_allRuns;
            varargout{4} = S_c_end_allRuns;
            if P.saveParSweepResToFile                        
                % Define name and save
                resName_base = ['./outputs/parameter_sweep_ASM1seq/results_t',num2str(P.T_sim),'_',char(datetime("today"))];
                filename_pSweepASM1seq_results = resName_base;
                parSetsName_base = ['./outputs/parameter_sweep_ASM1seq/parSets_t',num2str(P.T_sim),'_',char(datetime("today"))];
                parSetsName = parSetsName_base;
                cc = 1;
                while exist([filename_pSweepASM1seq_results,'.mat'],'file')
                    filename_pSweepASM1seq_results = [resName_base,'_v',num2str(cc)];
                    parSetsName = [parSetsName_base,'_v',num2str(cc)];
                    cc = cc+1;
                end
                save([filename_pSweepASM1seq_results,'.mat'],'results','grNet_AUT_allRuns','X_c_end_allRuns','S_c_end_allRuns','-v7.3')
                paramSetsUsed = P.pSweep.paramSetsUsed;
                save([parSetsName,'.mat'],'paramSetsUsed')
            end
        case 'EEMorris'
            fetchedRes = fetchOutputs(f); % If it does not work try to add 'UniformOutput',false
            sum(isnan(fetchedRes));
            results = reshape(fetchedRes,P.EE.k+1,P.EE.nTraj);  
            varargout{1} = results;
            save('outputs/parameter_Morris/morrisRuns_fetchedReshaped.mat','results');
            outputTable = array2table([usedParams;fetchedRes'],'VariableNames',cellstr(strcat('r',string(ii_all))));    
            writetable(outputTable,P.outputMorrisFileName,'Sheet','runs','Range','G1')
    end
 end
return

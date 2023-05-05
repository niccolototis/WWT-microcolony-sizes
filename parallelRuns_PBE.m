function [] = parallelRuns_PBE(aSt,P,M_Nti_0_full,t_all_ASM1seq,S_c_V_whole_tall,Ker_B,Ker_V,sp)
 
 delete(gcp('nocreate'))
 
 % - setups
 P.nSweeps = length(P.lambdasToSweep);
 P.figure = 'flat';
 P = defineColorMap(P);
 nRuns = numel(P.lambdasToSweep); 

    
 % Create a figure object, and set 'Visible' to true so that it opens in a new window, outside of the live script. To visualize the results of the parameter sweep, create a thisFig plot. Note that initializing the Z component of the thisFig with NaN creates an empty plot.
 switch P.figure
    %         case 'surface' 
    %             thisFig = surf(P.p1,P.p2,NaN(size(P.p1)));
    %             % note that the matrix for the Z coord needs to be given with dimensions in the opposite order
    %             % e.g: surf(1:5,1:2,rand(2,5))
    %             xlabel('paramLevel')
    %             ylabel('trajectory')
    %             zlabel('relative AUT growth')
    %             P.title = 'Morris (k+1) runs for each trajectory';
    %             title(P.title)                    
    case 'flat'
        thisFig = figure('Visible',true);
        P.lgdEntries = [];
        P.lgdfullNames = [];
        hold on
 end

 % Perform Parallel Parameter Sweep
 % refer to https://nl.mathworks.com/help/parallel-computing/run-code-on-parallel-pools.html#:~:text = What%20Is%20a%20Parallel%20Pool,cluster%20in%20your%20parallel%20preferences.
 % After you define the parameters, you can perform the parallel parameter sweep.
 % parfeval works more efficiently when you distribute the workload. To distribute the workload, group the parameters to explore into partitions. For this example, split into uniform partitions of size step by using the colon operator (:). The resulting array partitions contains the boundaries of the partitions. Note that you must add the end point of the last partition.
 tryThis = 1; % or 2
 stepPart = max(1,round(10^(numel(num2str(nRuns))-tryThis)));
 partitions = [1:stepPart:nRuns, nRuns+1]; % numel is the number of array elements
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
 parpool;

 % To send data from the workers, create a DataQueue object. Set up a function that updates the thisFig plot each time a worker sends data by using the afterEach function. The updateASM1seqPlot function is a supporting function defined at the end of the example.
 Q = parallel.pool.DataQueue;
 
 afterEach(Q,@(data) updateASM1seqPlotPBE(thisFig,data,P));
 

 for np = 1:numel(partitions)-1 % np=number of partition
    first = partitions(np);
    last = partitions(np+1);
    f(np,:) = parfeval(@collectResultOnePartitionPBE,1,first,last,aSt,P,M_Nti_0_full,t_all_ASM1seq,S_c_V_whole_tall,Ker_B,Ker_V,sp,Q);
 end

 % parfeval does not block MATLAB, so you can continue working while computations take place. The workers compute in parallel and send intermediate results through the DataQueue as soon as they become available.
 % If you want to block MATLAB until parfeval completes, use the wait function on the future objects. Using the wait function is useful when subsequent code depends on the completion of parfeval.
 wait(f);   
 % filename = strcat(['outputs/data/EE_parfeval_output_M',num2str(P.EE.M),'r',num2str(P.EE.r),'p',num2str(P.EE.p),'k',num2str(P.EE.k),'.mat']);

 objFns = fetchOutputs(f); % If it does not work try to add 'UniformOutput',false
 outputMatrix = [P.lambdasToSweep(:), objFns];
 figure
 plot(P.lambdasToSweep,objFns)
 objNames = {'objVal-notWeighted-N','objVal-notWeighted-pdf','objVal-notWeighted-pvf','objVal-weighted-N','objVal-weighted-pdf','objVal-weighted-pvf'};
 varNames = ['lambda' objNames];
 outputTable = array2table(outputMatrix,'VariableNames',varNames);
 dir_sweepLambda = 'outputs/parameter_sweep_PBE_lambda/';
 if ~exist(dir_sweepLambda)
    mkdir(dir_sweepLambda);
 end
 filename_sweepLambda = 'SweepLambda_objValues.xlsx';
 writetable(outputTable,[dir_sweepLambda filename_sweepLambda],'WriteVariableNames',true,'WriteMode','overwritesheet');                    
 spl = stackedplot(outputTable,objNames);
 spl.LineWidth = 2;
 spl.xlabel('Time [days]');
 spl.FontSize = 15;
 delete(gcp('nocreate'))
return




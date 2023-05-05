function [P,AllSweep_ASM1seq] = plotOneBatch_pSweep_ASM1seq(P,AllSweep_ASM1seq,qStr,X_HET_0_lb,X_HET_0_ub,X_AUT_0_lb,X_AUT_0_ub,X_c_res_base,X_c_res)
 results = AllSweep_ASM1seq.(['batch_',qStr]).results;
 nRuns = size(results,1);   
 runIdx = AllSweep_ASM1seq.(['batch_',qStr]).runIdx;
 X_c_end_allRuns = AllSweep_ASM1seq.(['batch_',qStr]).X_c_end_allRuns;
 S_c_end_allRuns = AllSweep_ASM1seq.(['batch_',qStr]).S_c_end_allRuns;
 
 % Runs_acceptance 5 columns
 % 1 - original runIdx
 % 2 - within allowd bounds for HET = 1, outside = 0
 % 3 - within allowd bounds for AUT = 1, outside = 0
 % 4 - computed metric on the distance between steady state vector and
 % initial condition vector
 % 5 - from the steadystate active biomass concentrations I compute the
 % total MLSS -> will be compared to the experimental data
 % 6 - ranking, 1 = run with lowest computed metric
 
 Runs_acceptance = array2table(nan(nRuns,6),'VariableNames',{'runIdx','inBounds_HET','inBounds_AUT','distance_SS_to_IS','total_MLSS_inBounds','Ranking'});
 Runs_acceptance.runIdx = runIdx;
 
 plotname = {'HET','AUT'};
 
 % save('./outputs/tmpOutputs/runASM1seqparamSweep_tmp.mat')

 % load('./outputs/tmpOutputs/runASM1seqparamSweep_tmp.mat')

 %%
 
 figHET = figure('visible','on','Position',[892,385,560,420]);
 figHET.Name = 'figHET';
 figAUT = figure('visible','on','Position',[881,54,560,420]);
 figAUT.Name = 'figAUT';

 % Removing NaNs
 NaN_rows = find(all(isnan(results),2));
 results(NaN_rows,:) = [];
 X_c_end_allRuns(NaN_rows,:) = [];
 S_c_end_allRuns(NaN_rows,:) = [];
 runIdx(NaN_rows,:) = [];
 Runs_acceptance(NaN_rows,:) = [];

 % Resizing the color map based on only the effective runs. In this way
 % I avoid, if I have only few runs acceptable, that the colors are too
 % similar
 nRuns = size(results,1);   
 P.nSweeps = nRuns;
 P = defineColorMap(P);

 for aRun = 1:nRuns
    lineHET = results(aRun,1:P.ntpts_ASM1seq);
    lineAUT = results(aRun,P.ntpts_ASM1seq+1:P.ntpts_ASM1seq*2);
    P = plotOneRunASM1seq(P,figHET,lineHET,'HET',[],'paramSweepHET',aRun);        
    P = plotOneRunASM1seq(P,figAUT,lineAUT,'AUT',[],'paramSweepAUT',aRun);

    %% computing acceptance score:
    % 1) It checks which of the runs have HET and AUT biomass that ends up in
    % the admissible region
    % 2) Checks that when I have 97% HET I do not have more than 3% AUT
    % 3) Computes the total relative distance
    % sum_i([(x_dev_i-x_ref_i)/x_ref_i]) between all the other steady state
    % vector variables (other than X_HET and X_AUT) and their values as
    % defined in defineInitialConditions. The reason for 3) is that I might
    % obtain a steady state vector at t_end of my simulation in which for
    % instance the value of nitrogen concentration is very far from the
    % initial one, so it might not be in line with the concentrations
    % observed experimentally at the steady state       

    X_c_end = X_c_end_allRuns(aRun,:)';
    S_c_end = S_c_end_allRuns(aRun,:)';
    if ~reallyEqual(X_c_end(P.HET),lineHET(end)) || ~reallyEqual(X_c_end(P.AUT),lineAUT(end))
        disp([X_c_end(P.HET),lineHET(end);X_c_end(P.AUT),lineAUT(end)])
        error('these should be the same')
    end

    % Check that HET_end is within bounds
    if X_c_end(P.HET) >= X_HET_0_lb && X_c_end(P.HET) <= X_HET_0_ub
        Runs_acceptance.inBounds_HET(aRun) = 1; % Within bounds = success!
    else
        Runs_acceptance.inBounds_HET(aRun) = 0; % Out of bounds
    end

    % Check that AUT_end is within bounds
    if X_c_end(P.AUT) >= X_AUT_0_lb && X_c_end(P.AUT) <= X_AUT_0_ub
        Runs_acceptance.inBounds_AUT(aRun) = 1; % Within bounds = success!
    else
        Runs_acceptance.inBounds_AUT(aRun) = 0; % Out of bounds
    end

    P_tmp = defineInletConcentrations(P);
    P_tmp = defineInitialConditions_t0(P_tmp);
    X_c_0 = P_tmp.this_x_all_0(1:P.nPart);
    S_c_0 = P_tmp.this_x_all_0(P.nPart+1:P.nPart+P.nSolub);
    [dist_ss_from_IC] = compare_and_sum(X_c_end,S_c_end,X_c_0,S_c_0);
    Runs_acceptance.distance_SS_to_IS(aRun) = dist_ss_from_IC;    
 
    X_active = X_c_end(P.HET)+X_c_end(P.AUT);
    X_MLSS_lb = X_active/(P.pSweep.tableBounds{'active_biom_in_gCOD','ub'}*P.pSweep.tableBounds{'gCOD_in_gVSS','ub'}*P.pSweep.tableBounds{'gVSS_in_gMLSS','ub'});
    X_MLSS_ub = X_active/(P.pSweep.tableBounds{'active_biom_in_gCOD','lb'}*P.pSweep.tableBounds{'gCOD_in_gVSS','lb'}*P.pSweep.tableBounds{'gVSS_in_gMLSS','lb'});

    % check if 7000 is within this bounds == can be explained with a
    % allowed combo of the other parameters here at the denominator
    Runs_acceptance.total_MLSS_inBounds(aRun) = (P.pSweep.ref_biom_conc>X_MLSS_lb & P.pSweep.ref_biom_conc<X_MLSS_ub);
 end
 % Create the rankings: Put all zeroes for the runs which fell out of
 % the bounds
 whereInBounds = (Runs_acceptance.inBounds_HET==1 & Runs_acceptance.inBounds_AUT==1 & Runs_acceptance.total_MLSS_inBounds==1);
 Runs_acceptance.Ranking(~whereInBounds) = 0; 
 % Sort the other ones based on distance from initial state (min to max)
 [~,idx] = sort(Runs_acceptance.distance_SS_to_IS(whereInBounds));
 ranking = zeros(size(idx));
 ranking(idx) = (1:length(idx));
 Runs_acceptance.Ranking(whereInBounds) = ranking;

 % Add the statistics on the 'goodness of each parameter set' to the
 % matrix

 %% Plotting baseline and chosen trajectory + an alternative
 set(figHET,'visible','on')
 set(figAUT,'visible','on')

 % - for baseline params
 P = plotOneRunASM1seq(P,figHET,X_c_res_base(:,P.HET),'HET','Nominal - literature','Nominal - literature',[]);
 P = plotOneRunASM1seq(P,figAUT,X_c_res_base(:,P.AUT),'AUT','Nominal - literature','Nominal - literature',[]);

 % - for chosen params
 P = plotOneRunASM1seq(P,figHET,X_c_res(:,P.HET),'HET','Chosen set','Chosen set',[]);
 P = plotOneRunASM1seq(P,figAUT,X_c_res(:,P.AUT),'AUT','Chosen set','Chosen set',[]);

% - Alternative set to the one chosen
aRun = find(Runs_acceptance{:,end}==1);
if ~isempty(aRun)
    oneRound = 1;
    oneLine = results(aRun,((oneRound-1)*P.ntpts_ASM1seq+1):(oneRound*P.ntpts_ASM1seq));
    P = plotOneRunASM1seq(P,figHET,oneLine,char(plotname(oneRound)),'Alternative set',[],aRun,'green');
    oneRound = 2;
    oneLine = results(aRun,((oneRound-1)*P.ntpts_ASM1seq+1):(oneRound*P.ntpts_ASM1seq));
    P = plotOneRunASM1seq(P,figAUT,oneLine,char(plotname(oneRound)),'Alternative set',[],aRun,'green');
end

 %% Plotting bands for admissible region
figure(figHET)
axes = gca;
hold on
r1 = fill(axes,[axes.XLim(1) axes.XLim(1) axes.XLim(end) axes.XLim(end)],[0 X_HET_0_lb X_HET_0_lb 0],'k','FaceAlpha',0.2,'EdgeColor','none'); % Plots only the shapes with no stroke   
set(get(get(r1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
uistack(r1,'bottom');
maxY = axes.YLim(end);
maxY = max(maxY,X_HET_0_ub+X_HET_0_ub/7);
r2 = fill(axes,[axes.XLim(1) axes.XLim(1) axes.XLim(end) axes.XLim(end)],[X_HET_0_ub maxY maxY X_HET_0_ub],'k','FaceAlpha',0.2,'EdgeColor','none'); % Plots only the shapes with no stroke   
set(get(get(r2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
uistack(r2,'bottom');
axes.YLim(end) = maxY;
hold off

figure(figAUT)
hold on
axes = gca;
r1 = fill(axes,[axes.XLim(1) axes.XLim(1) axes.XLim(end) axes.XLim(end)],[0 X_AUT_0_lb X_AUT_0_lb 0],'k','FaceAlpha',0.2,'EdgeColor','none'); % Plots only the shapes with no stroke   
set(get(get(r1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
uistack(r1,'bottom');
maxY = axes.YLim(end);
maxY = max(maxY,X_AUT_0_ub+X_AUT_0_ub/7);
r2 = fill(axes,[axes.XLim(1) axes.XLim(1) axes.XLim(end) axes.XLim(end)],[X_AUT_0_ub maxY maxY X_AUT_0_ub],'k','FaceAlpha',0.2,'EdgeColor','none'); % Plots only the shapes with no stroke   
set(get(get(r2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
uistack(r2,'bottom');
axes.YLim(end) = maxY;
hold off

 %% Adjust HET fig
 figure(figHET)
 axes = gca;
 axes.TickDir = 'out';
 axes.YTick = [0 1000 2000 3000 4000];           
 axes.YAxis.Exponent = 3;
 set(axes,'OuterPosition',[0.1 0.1 0.9 0.9]);
 axes.XLabel.Units = 'pixels';
 axes.YLabel.Units = 'pixels';
 axes.XLabel.FontSize = 20;
 axes.YLabel.FontSize = 20;
 axes.XLabel.Position(2) = axes.XLabel.Position(2)-5;
 axes.YLabel.Position(1) = axes.YLabel.Position(1)-10;
 % - Title
 axes.Title.FontSize = 25;
 axes.Title.Units = 'pixels';
 axes.Title.Position(2) = axes.Title.Position(2)+8;
 axes.Position(1) = 0.15;
 axes.Position(2) = 0.15;
 thisLeg = findobj(gcf, 'Type', 'Legend');
 thisLeg.FontSize = 18;

 %% Adjust AUT fig
 figure(figAUT)
 axes = gca;
 axes.TickDir = 'out';
 set(axes,'OuterPosition',[0.1 0.1 0.9 0.9]);
 axes.XLabel.Units = 'pixels';
 axes.YLabel.Units = 'pixels';
 axes.XLabel.FontSize = 20;
 axes.YLabel.FontSize = 20;
 axes.XLabel.Position(2) = axes.XLabel.Position(2)-5;
 axes.YLabel.Position(1) = axes.YLabel.Position(1)-10;
 % - Title
 axes.Title.FontSize = 25;
 axes.Title.Units = 'pixels';
 axes.Title.Position(2) = axes.Title.Position(2)+8;
 axes.Position(1) = 0.15;
 axes.Position(2) = 0.15;
 thisLeg = findobj(gcf, 'Type', 'Legend');
 thisLeg.FontSize = 18;

 %% SAVE
 orient(figHET,'landscape')
 saveas(figHET,['./outputs/figures/ASM1seq/parSweepASM1seq_',['batch_',qStr],'HET.pdf']);
 saveas(figHET,['./outputs/figures/ASM1seq/parSweepASM1seq_',['batch_',qStr],'HET.png']);
 
 orient(figAUT,'landscape')
 saveas(figAUT,['./outputs/figures/ASM1seq/parSweepASM1seq_',['batch_',qStr],'AUT.pdf']);
 saveas(figAUT,['./outputs/figures/ASM1seq/parSweepASM1seq_',['batch_',qStr],'AUT.png']);
 
 AllSweep_ASM1seq.X_c_end_allRuns = [AllSweep_ASM1seq.X_c_end_allRuns;X_c_end_allRuns];
 AllSweep_ASM1seq.S_c_end_allRuns = [AllSweep_ASM1seq.S_c_end_allRuns;S_c_end_allRuns];
 AllSweep_ASM1seq.Runs_acceptance = [AllSweep_ASM1seq.Runs_acceptance;Runs_acceptance];
 AllSweep_ASM1seq.runIdx = [AllSweep_ASM1seq.runIdx;runIdx];
return

function [] = plotLambdaProfile(statsToCompute,scatterWhatprofile,lambdaVs,howDistancesComputed,divideDistByCritValue,logDistance)
global lambdasToSweep folderOutput nDecimalsToIncludeNDFsampling 
% In the following code the commented for loops might be used in order to
% compute the distances among the simulations

 aStat = 'CVM';
 aW = 'equal';
 % Specific variables that are needed to
 nStat = find(strcmp(statsToCompute,scatterWhatprofile.stat));
 aW = scatterWhatprofile.weight;
 [foundLambdas,foundDists] = deal(nan(length(lambdasToSweep),1));
 switch lambdaVs
    case 'sum_dist_all'
        addToStat = {'','_all','_all'};
    case 'dist_tend'
        addToStat = {'','',''};
 end
 weights = {'equal','NDF','PDF'};
 r = 0;
 thisFig = figure('Position',[7,477,678,305]);
 hold on
 %    for nStat = 1:length(statsToCompute)
 %        aStat = statsToCompute{nStat};
 %        for nW = 1:length(weights)
 %            aW = weights{nW};
 %            optLam.(aStat).(aW).dist = Inf;
 %            r = r+1;
 %            thisFig = subscatter(length(statsToCompute),length(weights),r);
 %            hold on
 %            title([aStat ' ' aW])
 %        end
 %    end
 switch howDistancesComputed
    case 'tptsApart'
        load([folderOutput 'Dist_lambdasToSweep_Apart_' num2str(nDecimalsToIncludeNDFsampling) 'dec.mat'],'Dist_lambdasToSweep')
    case 'tptsAllin'     
        load([folderOutput 'Dist_lambdasToSweep_Allin_' num2str(nDecimalsToIncludeNDFsampling) 'dec.mat'],'Dist_lambdasToSweep')
 end
 % This just to load back P.lambdasToSweep   
 %  colori = {'r','y','g'};
 for i = 1:length(lambdasToSweep)   
    Lambda_i = lambdasToSweep(i);
    labLambda = ['Lambda_' num2str(Lambda_i)]; 
    labLambda = strrep(labLambda,'.','_');
%     for nStat = 1:length(statsToCompute)
        aStat = statsToCompute{nStat};
        aTS = addToStat{nStat};
        fprintf('Lambda %g p-vals for %s: ',Lambda_i,aStat);
%         for nW = 1:length(weights)
%             aW = weights{nW};
                try
                    if strcmp(lambdaVs,'sum_dist_all') || strcmp(aStat,'KS')
                        aDist = Dist_lambdasToSweep.(labLambda).(aW).mydist.([aStat aTS]);
                    else
                        aDist = Dist_lambdasToSweep.(labLambda).(aW).mydist.(aStat)(end);
                    end
                    if divideDistByCritValue && ~strcmp(aStat,'KS')
                        aDist = aDist/Dist_lambdasToSweep.(labLambda).(aW).crit.(aStat);
                    end
                    %                     if aDist<optLam.(aStat).(aW).dist
                    %                         optLam.(aStat).(aW).dist = aDist;
                    %                         optLam.(aStat).(aW).lambda = Lambda_i;
                    %                         optLam.(aStat).(aW).i = i;
                    %                     end
                    fprintf('%s = %g  ',aW,aDist);
                    foundDists(i) = aDist;
                    foundLambdas(i) = Lambda_i;
                    % scatter(thisFig,Lambda_i,aDist,'bo','MarkerFaceColor',colori{nStat});
                    % xline(thisFig,Lambda_i)
        %             scatter(Lambda_i,aDist,'bo','MarkerFaceColor',colori{nStat});
        %             xline(Lambda_i)
                catch
                end
        %end
        fprintf('\n');
    %end
 end
 
 fprintf('Adding missing lambdas \n');
 presentLambdas_idx = find(~isnan(foundLambdas));
 missingLambdas_idx = find(isnan(foundLambdas)); bb=find(isnan(foundDists));
 if ~isequal(missingLambdas_idx,bb)
    error('Should be the same')
 end
 missingLambdas_idx = sort(missingLambdas_idx);
 foundLambdas(missingLambdas_idx) = lambdasToSweep(missingLambdas_idx);
 x = foundLambdas(presentLambdas_idx);
 y = foundDists(presentLambdas_idx);
 xnew = foundLambdas(missingLambdas_idx);
 foundDists(missingLambdas_idx) = interp1(x, y, xnew);
 optLam.(aStat).(aW).dist = Inf;
 for i = 1:length(foundLambdas)
    if foundDists(i)<optLam.(aStat).(aW).dist
        optLam.(aStat).(aW).dist = foundDists(i);
        optLam.(aStat).(aW).lambda = foundLambdas(i);
        optLam.(aStat).(aW).i = i;
    end
    aDist = foundDists(i);
    if logDistance(nStat)
       aDist = log(aDist);
    end
    plot(foundLambdas(i),aDist,'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','none','MarkerSize', 7);
    rr = xline(foundLambdas(i));
    uistack(rr,'bottom');
 end


 % Plotting with a different color the point with the min distance
 %    for nStat = 1:length(statsToCompute)
 %        aStat = statsToCompute{nStat};
 %        for nW = 1:length(weights)
 %            aW = weights{nW};
                 aDist = optLam.(aStat).(aW).dist;
                 if logDistance(nStat)
                   aDist = log(aDist);
                 end
                 gg = plot(optLam.(aStat).(aW).lambda,aDist,'ro','MarkerFaceColor','r','MarkerSize', 7);
                 uistack(gg,'top');
 %        end
 %    end
 
 xlabel('$\lambda$','Interpreter','Latex','fontweight','bold')
 ylabel('$log(\bar{d})$','Interpreter','Latex','fontweight','bold')
 title('Mean Cramer-von Mises distance $\bar{d}$','Interpreter','Latex')
 theseaxes = gca();
 theseaxes.YLabel.Position(1) = theseaxes.YLabel.Position(1)-0.9;
 theseaxes.Title.Position(2) = theseaxes.Title.Position(2)+0.1;
 a = get(theseaxes,'XTickLabel');  
 set(theseaxes,'XTickLabel',a,'fontsize',15);
 a = get(theseaxes,'YTickLabel');  
 set(theseaxes,'YTickLabel',a,'fontsize',15);
 theseaxes.XLabel.FontSize = 20;
 theseaxes.YLabel.FontSize = 20;
 theseaxes.Title.FontSize = 20;
 theseaxes.OuterPosition = [0.01,0.01,1,1];

 foldFig = './outputs/figures/lambda_profile/';
 saveas(gcf,[foldFig 'lambda_profile_OK_',aStat,'_',lambdaVs,'.fig'])
 saveas(gcf,[foldFig 'lambda_profile_OK_',aStat,'_',lambdaVs,'.pdf'])
return

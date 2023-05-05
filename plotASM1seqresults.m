function [] = plotASM1seqresults(P,t_all,X_c_tall,S_c_tall,V_tall)
global gcs T_sim nnsos 
 if strcmp(num2str(T_sim),num2str(P.CycleLength))
    cc = '_1C';
 else
    cc = '';
 end
 ltys = {'-','o'};
 
 % Plot HET 
 colors = distinguishable_colors(P.nSp);
 figure
 for nnso = nnsos 
    P.gc = char(gcs(nnso));  
    lty = char(ltys(nnso));
    plot(t_all(1:100:end), X_c_tall.(P.gc)((1:100:end),P.HET),lty, 'LineWidth', 2,'Color',colors(P.HET,:))
    hold on 
 end
 for kk = 1:P.nPhases_1C
    xline(P.switchTimes_all(kk),'r-')
    hold on
 end
 hold off
 title('Biomass HET')
 lgdNames = P.spNames;
 lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
 lgdNames = lgdNames(:)';
 legend(lgdNames ,'fontsize',10);
 set(gca,'fontsize',15)
 xlabel('Time (days)') 
 ylabel('Concentration (mg/L) ')
 filename = strcat([P.folderPlots,'Biom',cc,'.pdf']);
 saveas(gcf,filename)
 
 % Plot AUT
 colors = distinguishable_colors(P.nSp);
 figure
 for nnso = nnsos 
    P.gc = char(gcs(nnso));  
    lty = char(ltys(nnso));
    plot(t_all(1:100:end), X_c_tall.(P.gc)((1:100:end),P.AUT),lty, 'LineWidth', 2,'Color',colors(P.AUT,:))
    hold on 
 end
 for kk = 1:P.nPhases_1C
    xline(P.switchTimes_all(kk),'r-')
    hold on
 end
 hold off
 title('Biomass AUT')
 lgdNames = P.spNames;
 lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
 lgdNames = lgdNames(:)';
 legend(lgdNames ,'fontsize',10);
 set(gca,'fontsize',15)
 xlabel('Time (days)') 
 ylabel('Concentration (mg/L) ')
 filename = strcat([P.folderPlots,'Biom',cc,'.pdf']);
 saveas(gcf,filename)
    
 % Plot particulate components 
 figure
 colors = distinguishable_colors(P.nPart+P.nSp);
 for aPart = (P.nSp+1):P.nPart
    for nnso = nnsos              
        P.gc = char(gcs(nnso));  
        lty = char(ltys(nnso));
        plot(t_all(1:100:end), X_c_tall.(P.gc)((1:100:end),aPart),lty,  'LineWidth', 2,'Color',colors(aPart,:))
        hold on
    end
 end
 for kk = 1:P.nPhases_1C
    xline(P.switchTimes_all(kk),'r-')
    hold on
 end
 hold off
 title('Particulate components')
 lgdNames = P.partNames;
 lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
 lgdNames = lgdNames(:)';
 legend(lgdNames ,'fontsize',10);
 set(gca,'fontsize',15)
 xlabel('Time (days)') 
 ylabel('Concentration (mg/L) ')
 hold off
 filename = strcat([P.folderPlots,'Particles',cc,'.pdf']);
 saveas(gcf,filename)
 
 % Plot soluble components  
 figure
 colors = distinguishable_colors(P.nSolub);
 for aSub = 1:P.nSolub
    for nnso = nnsos              
        P.gc = char(gcs(nnso));  
        lty = char(ltys(nnso));
        plot(t_all(1:100:end), S_c_tall.(P.gc)((1:100:end),aSub),lty, 'LineWidth', 2,'Color',colors(aSub,:))
        hold on
    end
 end
 for kk = 1:P.nPhases_1C
    xline(P.switchTimes_all(kk),'r-')
    hold on
 end
 hold off
 title('Soluble components')
 lgdNames = P.solNames;
 lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
 lgdNames = lgdNames(:)';
 legend(lgdNames ,'fontsize',10);
 set(gca,'fontsize',15)
 xlabel('Time (days)') 
 ylabel('Concentration (mg/L) ')
 hold off
 filename = strcat([P.folderPlots,'Solubles',cc,'.pdf']);
 saveas(gcf,filename)
 
 % Plot NH components  
 figure
 for nnso = nnsos             
    P.gc = char(gcs(nnso));  
    lty = char(ltys(nnso));
    plot(t_all(1:100:end), S_c_tall.(P.gc)((1:100:end),P.sNH4),lty, 'LineWidth', 2,'Color',colors(P.sNH4,:))
    hold on
    plot(t_all(1:100:end), S_c_tall.(P.gc)((1:100:end),P.sNO),lty, 'LineWidth', 2,'Color',colors(P.sNO,:))
    hold on
    plot(t_all(1:100:end), S_c_tall.(P.gc)((1:100:end),P.sND),lty, 'LineWidth', 2,'Color',colors(P.sND,:))
    hold on
 end
 for kk = 1:P.nPhases_1C
    xline(P.switchTimes_all(kk),'r-')
    hold on
 end
 hold off
 title('Nitrogen')
 lgdNames = {'NH4','NO','ND'};
 lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
 lgdNames = lgdNames(:)';
 legend(lgdNames ,'fontsize',10);
 set(gca,'fontsize',15)
 xlabel('Time (days)') 
 ylabel('Concentration (mg/L) ')
 hold off
 filename = strcat([P.folderPlots,'Nitrogen',cc,'.pdf']);
 saveas(gcf,filename)
 
% % Plot Requirements
% figure
% otherNames = ["tCOD","tN"];
% colors = distinguishable_colors(length(otherNames));
% for aV = 1:length(otherNames)
%     for nnso = nnsos 
%         P.gc = char(gcs(nnso));  
%         lty = char(ltys(nnso));  
%         plot(t_all(1:100:end), others.(P.gc).(otherNames(aV))(1:100:end),lty, 'LineWidth', 2,'Color',colors(aV,:))
%         hold on
%     end
% end
% for kk = 1:P.nPhases_1C
%     xline(P.switchTimes_all(kk),'r-')
%     hold on
% end
% hold off
% title('Requirements')    
% lgdNames = otherNames;
% lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
% lgdNames = lgdNames(:)';
% legend(lgdNames ,'fontsize',10);
% set(gca,'fontsize',15)
% xlabel('Time (days)') 
% ylabel('Concentration (mg/L) ')
% hold off
% filename = strcat([P.folderPlots,'Requirements',cc,'.pdf']);
% saveas(gcf,filename)
 
 % Plot V   
 figure
 for nnso = nnsos 
    P.gc = char(gcs(nnso));  
    lty = char(ltys(nnso));         
    plot(t_all(1:100:end), V_tall.(P.gc)((1:100:end)),lty,'LineWidth', 2)
    hold on
 end
 title('Reactor Volume')
 set(gca,'fontsize',15)
 xlabel('Time (days)') 
 ylabel('Volume (m^3) ')
 hold off
 filename = strcat([P.folderPlots,'Volume',cc,'.pdf']);
 saveas(gcf,filename)
 
 if ~(t_all(end)<P.CycleLength)
    figure
    tspanOnecycle = 0:P.dtWholebag:P.CycleLength;
    endCycle = find(strcmp(string(num2cell(t_all)),num2str(P.CycleLength)));
    if length(endCycle)>1
        endCycle = min(endCycle);
    end
    for nnso = nnsos 
        P.gc = char(gcs(nnso));  
        lty = char(ltys(nnso)); 
        plot(t_all(1:endCycle), V_tall.(P.gc)(1:endCycle,end),lty ,'LineWidth', 2)
        hold on
    end
    for kk = 1:P.nPhases_1C
        xline(P.switchTimes_all(kk),'r-')
        hold on
    end
    hold off
    title('Vol One cycle')
    lgdNames = 'Reactor Volume';
    lgdNames = [strcat(lgdNames,'concs'); strcat(lgdNames,'grams')];
    lgdNames = lgdNames(:)';
    legend(lgdNames ,'fontsize',10);
    set(gca,'fontsize',15)
    xlabel('Time (mins)') 
    ylabel('Volume (m^3) ')
    filename = strcat([P.folderPlots,'VolOneCycle',cc,'.pdf']);
    saveas(gcf,filename)
 end
end

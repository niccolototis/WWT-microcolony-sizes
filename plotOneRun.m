function [] = plotOneRun(P,aFig,BiomassTrajectory,nameBiom,namePlot)
 figure(aFig)
 hold on
 plot(P.tspallin,BiomassTrajectory,'LineWidth', 1,'Color','blue'); % Redraw the line of baseline params if too many have superimposed and is no more visible
 legend('fontsize',10,'Location','northeast');
 arr = get(aFig,'Children');
 indx = find(strcmp(get(arr,'Tag'),'legend'));
 lgdfullNames = arr(indx).String;
 for kk = 1:P.nPhases_1C
    hold on
    xl = xline(P.switchTimes_all(kk),'r-');
    set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
 end
 legend(lgdfullNames);
 set(gca,'fontsize',15)
 title(nameBiom);
 xlabel('Time (d)') 
 ylabel('Concentration (mg/L) ')
 hold off
 set(gcf,'Visible','on')
 switch namePlot
    case 'baseline'
        aName = strcat([nameBiom,'_BASE']);
        filename = strcat([P.folderPlots,aName,'.fig']);
    case 'chosenSetASM1seq'
        aName = strcat([nameBiom,'_scan_',strjoin(P.parToSweep,'_')]);
        filename = strcat([P.folderPlots,aName,'.fig']);
 end
 savefig(filename)
return

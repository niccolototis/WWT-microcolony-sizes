function figs_ASM1seq = adjustPlotAndSave_lowSRT(figs_ASM1seq,P,figName)
global iT noTrials
aFigure = figure(figs_ASM1seq.(['fig',figName]));
% - Ticks and grid
axes = gca;
axes.TickDir = 'out';
switch P.parSetup
 case 'sweep'
    switch figName
        case 'HET'
            axes.YTick = [0 1000 2000 3000 4000];            
            axes.YAxis.Exponent = 3;
        case 'AUT'
            axes.YTick = [axes.YTick(1) 50 100];
    end
    grid on
    axes.XGrid = 'off';
 case 'chosenSetASM1seq'
    switch P.IC_at
        case 'SS_computed'
            switch figName
                case 'HET'
                    axes.Position = [0.117,0.2,0.6975,0.7335];
                    axes.InnerPosition = [0.117,0.2,0.6975,0.7335];
                    axes.OuterPosition = [0,0.101,0.9,0.9];
                    axes.YLim = [400 1900];
                    %
                    axes.YTick = [400:400:1600];      
                    axes.YTickLabel = axes.YTick/1000;
                    % axes.YAxis.Exponent = 3;
                case 'AUT'
                    axes.Position = [0.217,0.199,0.6975,0.7335];
                    axes.InnerPosition = [0.217,0.199,0.6975,0.7335];
                    axes.OuterPosition = [0.1,0.1,0.9,0.9];
                    axes.YLim = [18 74];
                    axes.YTick = [20 40 60];
            end
        case 'startupScenario'
            axes.Legend.Visible = 'off';
    end
end

if P.plotPhases
 %        for kk = 1:P.nPhases_1C
 %            hold on
 %            xl = xline(P.switchTimes_all(kk),'k-');
 %            set(get(get(xl,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
 %        end
 %P.phaseNames = {'onAer','','onFeed','onAer','onSettling','onPurge','onEfflux'}; 
 % P.phaseSeq_1C = [1 2 3 4 2 3 4 2 3 4 5 6 7]; % refers to the index of phaseNames
 axes = gca;

 % Here I am changing the limits because I assume that when I have
 % plotPhases I am plotting just the first cycle
 switch P.IC_at
    case 'SS_computed'
    case 'data'
        switch figName
            case 'HET'
                axes.YLim = [400 1900];
                axes.YTick = [ 1600 1800 2000 ];             
                axes.YAxis.Exponent = 3;
            case 'AUT'
                axes.YLim = [20 34];
                axes.YTick = [22 24 26 28 30 32];
        end
 end
 h = findobj(gca,'Type','line');
%     h.Color = [0 0 0];
 allObj = h;
 phaseNames = P.phaseNames(find(~cellfun(@isempty,P.phaseNames)));
 phaseNames = unique(phaseNames,'stable');
 phaseNamesLeg = strrep(phaseNames,'on','');
 phaseNamesLeg = strrep(phaseNamesLeg,'Aer','Aeration');
 nColors = 20;
 FaceColors = distinguishable_colors(nColors);
 fAlpha = 0.2;
 colorIdxForPhases = [10,2,18,9,3];
 %        % Try out the colors to choose them
 %        figure
 %        hold on
 %        for pp = 1:nColors
 %            b = bar(pp,10,'FaceColor',FaceColors(pp,:),'FaceAlpha',fAlpha);
 %        end
 %        %     {'Aeration'} - 10
 %        %     {'Feed'    } - 2
 %        %     {'Settling'} - 18
 %        %     {'Purge'   } - 9
 %        %     {'Efflux'  } - 3

 for zz = 1:length(phaseNames) % vado al contrario cosi viene giusta la legenda
    pos1 = find(strcmp(P.phaseNames,phaseNames{zz}));
    pos2 = find(ismember(P.phaseSeq_1C,pos1));
    axes = gca;
    yy = axes.YLim;
    Y = [yy(1) yy(2) yy(2) yy(1)];
    for gg = 1:length(pos2)
        X = [P.switchTimes_all(pos2(gg)) P.switchTimes_all(pos2(gg)) P.switchTimes_all(pos2(gg)+1) P.switchTimes_all(pos2(gg)+1)];
        ar.(phaseNames{zz}) = area(X,Y);
        uistack(ar.(phaseNames{zz}),'bottom') 
        ar.(phaseNames{zz}).LineStyle = 'none';
        ar.(phaseNames{zz}).FaceColor = FaceColors(colorIdxForPhases(zz),:);
        ar.(phaseNames{zz}).FaceAlpha = fAlpha;
        if gg ==1
            allObj = [allObj ar.(phaseNames{zz})];
            ar.(phaseNames{zz}).DisplayName = phaseNamesLeg{zz};
        else
            set(get(get(ar.(phaseNames{zz}),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
    end
 end
 switch P.IC_at
    case 'SS_computed'
        thisLeg = legend(allObj,'Location','southwest');
    case 'data'
        thisLeg = legend(allObj,'Location','southeast');
 end
end
thisLeg.FontSize = 18;

% - Axes labels position
f.Position = [20 50 600 400];
% f.OuterPosition(3) = f.OuterPosition(3)+100;
% f.OuterPosition(4) = f.OuterPosition(4)+50;

axes.FontSize = 32;

axes.XLim = [0 50];

set(axes,'OuterPosition',[0.1 0.1 0.9 0.9]);
axes.XLabel.Units = 'pixels';
axes.YLabel.Units = 'pixels';
axes.XLabel.FontSize = 32;
axes.YLabel.FontSize = 32;
axes.XLabel.Position(2) = axes.XLabel.Position(2)-2;
axes.YLabel.Position(1) = axes.YLabel.Position(1)-4;

% - Title
axes.Title.FontSize = 32;
axes.Title.Units = 'pixels';
axes.Title.Position(2) = 348;



hold off
set(aFigure,'Visible','on')

% - save fig
if iT ==noTrials
 P.allObjNames(end) = strrep(P.allObjNames(end),'d;','d.'); 
 
 [hleg,hobj] = legend(P.allObjNames,'interpreter','latex','location','northeast','Orientation','horizontal');   
 hleg.FontSize = 20;
 hleg.Position(2) = hleg.Position(2)-0.02;
 for p = (noTrials+1)+1:2:3*(noTrials+1)-1
    hobj(p).XData(1) = hobj(p).XData(1)+0.05;       
    hobj(p).YData = hobj(p).YData+0.1;
 end
 for p = (noTrials+1):-1:1
    hobj(p).FontSize = 28;  
    hobj(p).Units = 'normalized'; 
 end
 
 aa_ll = 0;
 for p = (noTrials+1):-1:1
    aa_ll = aa_ll+0.006;
    hobj(p).Position(1) = hobj(p).Position(1)+aa_ll;
 end
 aa_ll = 0;
 for p = 3*(noTrials+1)-1:-2:(noTrials+1)+1
    aa_ll = aa_ll+0.006;
    hobj(p).XData = hobj(p).XData+aa_ll;
    hobj(p).LineWidth = 20;
 end
 title(hleg,'$\mathrm{\mathbf{SRT = }}$','Interpreter','latex','FontSize',hobj(1).FontSize);
 hleg.Title.Visible = 'on';
 hleg.Title.NodeChildren.Position = [-0.08 0.52 0.6];
 legend('boxoff')
 
 set(gcf, 'Renderer', 'painters')

 saveNamePlot = [P.folderPlots,P.parSetup,'_',figName,'_lowSRT'];
 if strcmp(P.IC_at,'startupScenario')
    savefig([saveNamePlot ,P.IC_at,'.fig'])
    saveas(aFigure,[saveNamePlot,P.IC_at,'.pdf'])
 else
    saveas(gcf,['./outputs/figures/ASM1seq/',P.parSetup,'_',figName,'_lowSRT.pdf'])
    if P.plotPhases 
        savefig(aFigure,['./outputs/figures/ASM1seq/oneCycleSS_' figName '_lowSRT.fig'])
        saveas(aFigure,['./outputs/figures/ASM1seq/oneCycleSS_' figName '_lowSRT.pdf'])
    end
 end
end

figs_ASM1seq.(['fig' figName]) = aFigure;
end

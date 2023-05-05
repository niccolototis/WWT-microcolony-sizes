function varargout = adjustPlotAndSave(P,aFigure,figName,varargin)
global iT noTrials
figure(aFigure);
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
                    axes.YLim = [1460 1640];
                    axes.YTick = [1460 1500 1540 1580 1620];             
                    axes.YAxis.Exponent = 3;
                case 'AUT'
                    axes.YLim = [56 64];
                    axes.YTick = [56 58 60 62 64];
            end
        case 'data'
            switch figName
                case 'HET'
                    axes.YLim = [1400 2300];
                    axes.YTick = [1400 1600 1800 2000 2200];             
                    axes.YAxis.Exponent = 3;
                case 'AUT'
                    axes.YLim = [22 67];
                    axes.YTick = [20:10:70];
            end
        case 'startupScenario'
            axes.Legend.Visible = 'off';
            axes.XLim = [0 50];
            switch figName
                case 'HET'
                    axes.YTick = [ 500 1000 1500 ];             
                    axes.YAxis.Exponent = 3;
                case 'AUT'
                    axes.YLim = [0 70];
                    axes.YTick = [20 40 60];
            end
    end
end
axes.XLim = [0 P.T_sim];
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
                axes.YLim = [1500 2100];
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
set(axes,'OuterPosition',[0.1 0.1 0.9 0.9]);
axes.XLabel.Units = 'pixels';
axes.YLabel.Units = 'pixels';
axes.XLabel.FontSize = 20;
axes.YLabel.FontSize = 20;
axes.XLabel.Position(2) = axes.XLabel.Position(2)-2;
axes.YLabel.Position(1) = axes.YLabel.Position(1)-4;

% - Title
axes.Title.FontSize = 20;
axes.Title.Units = 'pixels';
axes.Title.Position(2) = axes.Title.Position(2)+10;


hold off
set(gcf,'Visible','on')

% if strcmp(figName,'HET')
%     thisLeg.Visible = 'off';
% end

% - save fig
saveNamePlot = [P.folderPlots,P.parSetup,'_',figName];
if nargin ==4
 saveNamePlot = [P.folderPlots,P.parSetup,'_',figName,varargin{1}];
end
if strcmp(P.IC_at,'startupScenario')
 axes.Title.FontSize = 32;
 axes.XLabel.FontSize = 32;
 axes.YLabel.FontSize = 32;
 axes.FontSize = 32;
 axes.XLabel.Position(2) = axes.XLabel.Position(2)-22;
 
 savefig([saveNamePlot ,P.IC_at,'.fig'])
 saveas(gcf,[saveNamePlot,P.IC_at,'.pdf'])
 saveas(gcf,['./outputs/figures/ASM1seq/',P.parSetup,'_',figName,'_',P.IC_at,'.pdf'])
else
 savefig([saveNamePlot,'.fig'])
 saveas(gcf,[saveNamePlot,'.pdf'])
 if P.plotPhases 
    savefig(gcf,['./outputs/figures/ASM1seq/oneCycleSS_' figName '.fig'])
    saveas(gcf,['./outputs/figures/ASM1seq/oneCycleSS_' figName '.pdf'])
 end
end

P.(['fig' figName]) = aFigure;
varargout{1} = P;
end

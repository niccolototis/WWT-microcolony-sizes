function [] = plotSSforPaper(P,BiomassTrajectories,namesBiom)
% - Ticks and grid
phaseNames = P.phaseNames(find(~cellfun(@isempty,P.phaseNames)));
phaseNames = unique(phaseNames,'stable');
phaseNamesLeg = strrep(phaseNames,'on','');
phaseNamesLeg = strrep(phaseNamesLeg,'Aer','Aeration');
nColors = 20;
FaceColors = P.allColors.phaseColors;
fAlpha = 0.3;

thisFig = figure('Position',[50,265,485,534]);
hold on
theseaxes = gca;
% - draw lines
idx1C = find(P.tspallin<=P.CycleLength);
dataHET = BiomassTrajectories(idx1C,P.HET);
lineHET = plot(P.tspallin(idx1C),dataHET,'LineWidth', 3,'Color',P.allColors.baseHET); % Redraw the line of baseline params if too many have superimposed and is no more visible
dataAUT = BiomassTrajectories(idx1C,P.AUT);
howManyRanges = 1.8;
dataAUTconverted = convertToFakeCoord(dataAUT,dataHET,dataAUT,howManyRanges);
lineAUT = plot(P.tspallin(idx1C),dataAUTconverted,'LineWidth', 3,'Color',P.allColors.baseAUT); % Redraw the line of baseline params if too many have superimposed and is no more visible
ticksHet = [1450 1500 1550 1600 1650]; 
ticksAUT = [ 56 58 60 62 64];
ticksAUT_converted = convertToFakeCoord(ticksAUT,dataHET,dataAUT,howManyRanges);
theseaxes.YTick = [ticksAUT_converted ticksHet];
theseaxes.YTickLabel = [string(ticksAUT) string(ticksHet)];
newmin = convertToFakeCoord(52,dataHET,dataAUT,howManyRanges);
oldLims = sort([newmin 1700]);
theseaxes.YLim = [min(oldLims) max(oldLims)];
x = yline(1430,'--k');
set(get(get(x,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

newcoor = convertToFakeCoord(65,dataHET,dataAUT,howManyRanges);
xx = yline(newcoor,'--k');
set(get(get(xx,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


%###################
% PLOT BANDS   %
%###################
%     %     {'Aeration'} - 10
%     %     {'Feed'    } - 2
%     %     {'Settling'} - 18
%     %     {'Purge'   } - 9
%     %     {'Efflux'  } - 3

theseaxes = gca();
%     h.Color = [0 0 0];
allObj = [lineHET lineAUT];
allObjNames = {'$x_{\mathrm{H}}$','$x_{\mathrm{A}}$'};
for zz = 1:length(phaseNames) % vado al contrario cosi viene giusta la legenda
 pos1 = find(strcmp(P.phaseNames,phaseNames{zz}));
 pos2 = find(ismember(P.phaseSeq_1C,pos1));
 for gg = 1:length(pos2)
    X = [P.switchTimes_all(pos2(gg)) P.switchTimes_all(pos2(gg)) P.switchTimes_all(pos2(gg)+1) P.switchTimes_all(pos2(gg)+1)];
    Y = [theseaxes.YLim(1) theseaxes.YLim(end) theseaxes.YLim(end) theseaxes.YLim(1)];
    ar.(phaseNames{zz}) = area(theseaxes,X,Y);
    uistack(ar.(phaseNames{zz}),'bottom') 
    ar.(phaseNames{zz}).LineStyle = 'none';
    ar.(phaseNames{zz}).FaceColor = FaceColors{zz};
    ar.(phaseNames{zz}).FaceAlpha = fAlpha;
    if gg ==1
        allObjNames = [allObjNames phaseNamesLeg{zz}];
        allObj = [allObj ar.(phaseNames{zz})];
    else
        set(get(get(ar.(phaseNames{zz}),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
 end
end
thisLegend = legend(allObj,allObjNames,'Interpreter','latex','location','southwest');


% - adjust aesthetics
xlabel(['$\mathrm{\mathbf{Time \; (d)}}$'],'interpreter','latex')
ylabel(['$\mathrm{\mathbf{Concentration \; (gm}}\mathbf{^{-3})}$'],'interpreter','latex') 
theseaxes.TickDir = 'out';


%% Adjust and save
theseaxes.XLabel.Units = 'pixels';
theseaxes.YLabel.Units = 'pixels';
theseaxes.FontSize = 15;
thisLegend.FontSize = 15;
theseaxes.XLabel.FontSize = 20;
theseaxes.YLabel.FontSize = 20;
theseaxes.XLabel.Position(2) = theseaxes.XLabel.Position(2)-2;
theseaxes.YLabel.Position(1) = theseaxes.YLabel.Position(1)-4;

% - Title
theseaxes.Title.FontSize = 20;
theseaxes.Title.Units = 'pixels';
theseaxes.Title.Position(2) = theseaxes.Title.Position(2)+10;

hold off
set(gcf,'Visible','on')

% - save fig
savefig(gcf,['./outputs/figures/ASM1seq/oneCycleSS.fig'])
saveas(gcf,['./outputs/figures/ASM1seq/oneCycleSS.pdf'])
return

function [dataAUTconverted] = convertToFakeCoord(dataAUTtoConvert,dataHET,dataAUT,howManyRanges)
rangeHET = range(dataHET); 
rangeAUT = range(dataAUT);
minAUT = min(dataAUT);
minAUT_new = min(dataHET)-howManyRanges*rangeHET;
dataAUTconverted = ((dataAUTtoConvert-minAUT)/rangeAUT)*rangeHET+minAUT_new;
return

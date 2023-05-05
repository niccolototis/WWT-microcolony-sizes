function figs_ASM1seq = plotOneRunASM1seq_lowSRT(figs_ASM1seq,P,BiomassTrajectory,nameBiom)
 global iT
 aFigure = figure(figs_ASM1seq.(['fig' nameBiom]));
 hold on
 theseAxes = aFigure.CurrentAxes;
 transparency = 2;
 if isempty(theseAxes)
    theseAxes = axes('Parent', aFigure);
 end

 lwidth = 2;
 switch iT
    case 0
        %FaceAlpha = 0.4;
        %col = cmap(thisRowCol,:);
        col = hex2rgb('85C43C');
    case 1
        %FaceAlpha = 0.6;
        %col = hex2rgb('FD6039');
        col = hex2rgb('C06A72');
    case 2
        %FaceAlpha = 0.8;
        col = hex2rgb('7888C0');
    case 3
        %FaceAlpha = 1;
        col = hex2rgb('C09954');
 end
 % switch nameBiom
 %    case 'HET'
 %        switch iT
 %            case 0
 %                col = hex2rgb('AAE1B0');
 %            case 1
 %                col = hex2rgb('80D289');
 %            case 2
 %                col = hex2rgb('55C361');
 %            case 3
 %                col = hex2rgb('2BB43A');
 %        end
 %    case 'AUT'
 %        switch iT
 %            case 0
 %                col = hex2rgb('F5CEC5');
 %            case 1
 %                col = hex2rgb('EFB6A8');
 %            case 2
 %                col = hex2rgb('EB9D8B');
 %            case 3
 %                col = hex2rgb('E5856E');
 %        end
 % end
 
 % - draw line
 hold(theseAxes,'on')
 if P.plotPhases
    idx1C = find(P.tspallin<=P.CycleLength);
    thisLine = plot(theseAxes,P.tspallin(idx1C),BiomassTrajectory(idx1C),'LineWidth', lwidth,'Color',col); % Redraw the line of baseline params if too many have superimposed and is no more visible
 else
    thisLine = plot(theseAxes,P.tspallin,BiomassTrajectory,'LineWidth', lwidth,'Color',col); % Redraw the line of baseline params if too many have superimposed and is no more visible
 end
 hold(theseAxes,'off')
 
% if ~isempty(nameLine) 
%     arr = get(aFigure,'Children');
%     indx = find(strcmp(get(arr,'Tag'),'legend'));
%     if isempty(indx)
%         lgdLines = thisLine;
%         lgdStrings = nameLine;        
%     else
%         lgdStrings = arr(indx).String;
%         lgdStrings{end} = nameLine;   % Needs to be a row
%         lgdLines = aFigure.Children(2).Children;
%         lgdLines = lgdLines(1:length(lgdStrings));
%         lgdLines = flip(lgdLines);
%     end
%     legend(lgdLines,lgdStrings,'fontsize',15,'Location','northeast');
% end
 
 %% - adjust aesthetics
 figure(aFigure)
 set(theseAxes,'fontsize',15)
 switch nameBiom
    case 'HET'
        title('$\mathbf{x_{\mathrm{\mathbf{H}}} \; (\mathrm{\mathbf{gm}}^{-3})}$','Interpreter','latex');
    case 'AUT'
        title('$\mathbf{x_{\mathrm{\mathbf{A}}} \; (\mathrm{\mathbf{gm}}^{-3})}$','Interpreter','latex');
 end
 xlabel(['$\mathrm{\mathbf{Time \; (d)}}$'],'interpreter','latex')
 theseAxes = gca;
 theseAxes.FontSize = 18;
 theseAxes.Title.FontSize = 20;
 
 figs_ASM1seq.(['fig' nameBiom]) = aFigure;
return

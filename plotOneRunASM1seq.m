function P = plotOneRunASM1seq(P,figHandle,BiomassTrajectory,nameBiom,nameLine,saveNamePlot,aRun,varargin)
 global iT
 theseAxes = figHandle.CurrentAxes;
 transparency = 2;
 if isempty(theseAxes)
    theseAxes = axes('Parent', figHandle);
 end
 if isempty(nameLine)
    nameLine = 'null';
 end
 switch nameLine
    case 'Nominal - literature'
        col = P.allColors.baseLiterature;
        lwidth = 2;
    case 'Chosen set'
        lwidth = 2;
        switch nameBiom
            case 'HET'
                col = P.allColors.baseHET;
            case 'AUT'
                col = P.allColors.baseAUT;
        end
    case 'Alternative set'
        col = varargin{1};
        lwidth = 2;
    otherwise
        transparency = 0.1;
        switch saveNamePlot
            case 'paramSweepHET'
                col = P.allColors.HET(aRun,:);
            case 'paramSweepAUT'
                col = P.allColors.AUT(aRun,:);
        end
        lwidth = 1.5;
        nameLine = [];
 end
 
 % - draw line
 hold(theseAxes,'on')
 if contains(saveNamePlot,'paramSweep')
    if P.plotPhases
        idx1C = find(P.tspallin<=P.CycleLength);
        thisLine = patchline(theseAxes,P.tspallin(idx1C),BiomassTrajectory(idx1C),'LineWidth', lwidth,'edgecolor',col,'edgealpha',transparency); % Redraw the line of baseline params if too many have superimposed and is no more visible
    else
        thisLine = patchline(theseAxes,P.tspallin,BiomassTrajectory,'LineWidth', lwidth,'edgecolor',col,'edgealpha',transparency); % Redraw the line of baseline params if too many have superimposed and is no more visible
    end
 else
    if P.plotPhases
        idx1C = find(P.tspallin<=P.CycleLength);
        thisLine = plot(theseAxes,P.tspallin(idx1C),BiomassTrajectory(idx1C),'LineWidth', lwidth,'Color',col); % Redraw the line of baseline params if too many have superimposed and is no more visible
    else
        thisLine = plot(theseAxes,P.tspallin,BiomassTrajectory,'LineWidth', lwidth,'Color',col); % Redraw the line of baseline params if too many have superimposed and is no more visible
    end
 end
 hold(theseAxes,'off')
 
 if ~isempty(nameLine) 
    arr = get(figHandle,'Children');
    indx = find(strcmp(get(arr,'Tag'),'legend'));
    if isempty(indx)
        lgdLines = thisLine;
        lgdStrings = nameLine;        
    else
        lgdStrings = arr(indx).String;
        lgdStrings{end} = nameLine;   % Needs to be a row
        lgdLines = figHandle.Children(2).Children;
        lgdLines = lgdLines(1:length(lgdStrings));
        lgdLines = flip(lgdLines);
    end
    legend(lgdLines,lgdStrings,'fontsize',15,'Location','northeast');
 end
 
 %% - adjust aesthetics
 figure(figHandle)
 set(theseAxes,'fontsize',15)
 switch nameBiom
    case 'HET'
        title(['$\mathbf{x_{\mathrm{\mathbf{H}}} \; (\mathrm{\mathbf{gm}}^{-3}}$)'],'interpreter','latex')
    case 'AUT'
        title(['$\mathbf{x_{\mathrm{\mathbf{A}}} \; (\mathrm{\mathbf{gm}}^{-3}}$)'],'interpreter','latex')
 end
 xlabel(['$\mathrm{\mathbf{Time \; (d)}}$'],'interpreter','latex') 
 theseAxes = gca;
 theseAxes.FontSize = 18;
 theseAxes.Title.FontSize = 20;
return

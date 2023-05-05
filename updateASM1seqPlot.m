function [] = updateASM1seqPlot(data,P)
% Define another helper function that updates the thisFig plot when new data arrives.
if ~isnan(data(2))
 switch P.parSetup
    case 'EEMorris'
        disp(data(2))
        P.figMorris.ZData(data(1)) = data(2);    
        drawnow('limitrate');
    case 'paramSweep'
        rounds = 1;
        plotname = cellstr(P.plot);
        if strcmp(P.plot,'both_HET_AUT')
            if strcmp(P.figure,'flat') 
                rounds = 2;
                plotname = [{'HET'},{'AUT'}];
                if length(data(2:end))~= P.ntpts_ASM1seq*2
                    ME = MException('I should have the trajectories for AUT and HET concatenated one below the other');
                    throw(ME);
                end
            else
                ME = MException('I cannot plot both_HET_AUT on the surface plot');
                throw(ME);    
            end
        end
        switch P.figure
            case 'surface'
                figname = strcat(['fig',char(plotname)]);
                thisFig = figure(P.(figname));
                disp(data(2))
                thisFig.ZData(data(1)) = data(2);    
                drawnow('limitrate');
            case 'flat'
                thisLev = data(1);
                for oneRound = 1:rounds
                    thisColor = P.allColors.(char(plotname(oneRound)))(thisLev,:);
                    figname = strcat(['fig',char(plotname(oneRound))]);
                    thisFig = figure(P.(figname));
                    hold on
                    skip = 1;
                    oneLine = data(((oneRound-1)*P.ntpts_ASM1seq+1+skip):(oneRound)*P.ntpts_ASM1seq+skip);
                    thisLine = plot(P.tspallin,oneLine,'LineWidth', 1,'Color',thisColor);
%                     set(get(get(thisLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%                     lgd = strcat(['r',num2str(thisLev)]);
%                     lgd = [];
%                     for jj = 1:length(P.pSweep.parNamesSweep) % for each of the parameters that I change
%                         E = strcat(P.pSweep.parNamesSweep{jj},'=');
%                         str = strcat(['p' num2str(jj)]);
%                         E = strcat([' ',E,sprintf('%.4g', P.(str)(thisLev))]); %I assign the ii-th par value to its specific parameter name 
%                         lgd = [lgd E];
%                     end
%                     arr = get(thisFig,'Children');
%                     indx = find(strcmp(get(arr,'Tag'),'legend'));
%                     if ~isempty(indx)
%                         lgdfullNames = arr(indx).String;
%                         if P.suppressLegend
%                             lgdfullNames = lgdfullNames(1);   % Needs to be a row
%                         else
%                             lgdfullNames(end) = cellstr(lgd);   % Needs to be a row
%                         end
%                         lgd = lgdfullNames;
%                     end
%                     legend(lgd ,'fontsize',10,'Location','westoutside');
%                     drawnow('limitrate');
%                     hold off
                end
        end
 end
end
return

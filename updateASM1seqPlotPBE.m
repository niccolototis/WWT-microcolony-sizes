function [] = updateASM1seqPlotPBE(thisFig,data,P)
    thisLev = data(1);
    fprintf('Plotting run number %g \n',thisLev);
    typeLine = {'mo','c+','r*','g.','bx','k*'};
    hold on
    for j = 2:length(data)
        plot(data(1),data(j),typeLine{j-1}); 
    end
    lgd = strcat(['r',num2str(thisLev)]);
    arr = get(thisFig,'Children');
    indx = find(strcmp(get(arr,'Tag'),'legend'));
    if ~isempty(indx)
        lgdfullNames = arr(indx).String;
        if P.suppressLegend
            lgdfullNames = lgdfullNames(1);   % Needs to be a row
        else
            lgdfullNames(end) = cellstr(lgd);   % Needs to be a row
        end
        lgd = lgdfullNames;
    end
    legend(lgd ,'fontsize',10,'Location','westoutside');
    drawnow('limitrate');
    return
end

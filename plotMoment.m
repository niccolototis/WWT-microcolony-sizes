function [] = plotMoment(tspan,moments,whatMoment,folderName)
global spNames speciesPBE
switch whatMoment
 case 1
    mo = 'First';
 case 2
    mo = 'Second';
 case 3
    mo = 'Third';
end
figure
for aSpPBE = speciesPBE
 aSpN = char(spNames{aSpPBE});
plot(tspan,moments.x2.norm.(mo).(aSpN),'LineWidth', 4)
hold on
end
plotname = strcat([mo,' moment m_',num2str(whatMoment)]);
title(plotname,'FontSize', 18);
lgd = legend(spNames(speciesPBE));
lgd.FontSize = 18;
xlabel('Time [-]','FontSize', 18) 
ylabel(strcat(['m_',num2str(whatMoment),' [-]']),'FontSize', 25)
cax = gca;
xtick = linspace(cax.XLim(1),cax.XLim(2),5);
ytick = linspace(cax.YLim(1),cax.YLim(2),5);
set(gca, 'XTick',xtick , 'XTickLabel',xtick,'FontSize',18)
set(gca, 'YTick',ytick , 'YTickLabel',ytick,'FontSize',18)
filename = strcat(['M',num2str(whatMoment)]);
saveas(gcf,strcat([folderName,filename,'.pdf']))
return

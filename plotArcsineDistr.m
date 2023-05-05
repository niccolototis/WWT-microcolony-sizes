close all
thisFig = figure('Position',[658,162,916,643]);
forPresentaz = false;
xSc = 0.05;
shift = xSc;
x = sort(0:1/500:1);
y = 1./(pi * sqrt(x.*(1-x)));
x = x+shift;
hold on
% plot(x,y,'r-','linewidth',1)
aa = area(x,y,'LineStyle','none');
if ~forPresentaz
 aa.FaceColor = hex2rgb('E5856E');
else
 aa.FaceColor = [101 170 113]/255;
end
aa.EdgeColor = 'w';

xMinChild = xSc;
xMaxChild = xSc+1;
xMom = xSc+1+xSc;
mom = xline(xMom,'--');

rect1 = patch([0 0 xMinChild xMinChild],[0,10,10,0],[200 200 200]/255,'LineStyle','none');
rect1.FaceColor = 'k';
rect1.FaceAlpha = 0.1;
rect1 = candystripe(rect1,'Color','w','Width',2);

rect2 = patch([xMaxChild xMaxChild xMom xMom],[0,10,10,0],[200 200 200]/255,'LineStyle','none');
rect2.FaceColor = 'k';
rect2.FaceAlpha = 0.1;
rect2 = candystripe(rect2,'Color','w','Width',2);

expSc = log(xMinChild);
expMom = log(xMom);
expTicks = linspace(expSc,expMom,6);
diffz = unique(diff(expTicks));
expTicks = [expTicks expTicks(end)+diffz]; 
ticks = exp(expTicks);
ticks = [0 ticks];
for jj = 2:length(ticks)
 xline(ticks(jj),'--','linewidth',1)
end
% minChild = xline(xMinChild,'k','linewidth',1);
% maxChild = xline(xMaxChild,'k','linewidth',1);

xlim([0 xSc+1+xSc+0.2])
ylim([0 7])
yticks([])
xticks(ticks)

dd = gca;
dd.TickDir = 'out';
dd.TickLabelInterpreter = 'latex';
dd.FontSize = 26;
dd.XTickLabel = ["$0$","$s_1$","$s_2$","$s_3$","$s_4$","$s_5$","$\zeta$"];
dd.XAxis.MajorTickChild.LineWidth = 1;



if forPresentaz

 dd.XAxis.LineWidth = 3;
 dd.YAxis.LineWidth = 3;

 dd.FontSize = 22;
 dd.TickLabelInterpreter = 'latex';

 dd.XTickLabel = ["$0$","$s_1$","","","$s$","","$\zeta$"];
 % dd.XTickLabel = cellstr(strcat("s_",string(1:length(ticks))));

 % for j = 1:length(ticks)
 %    plot([ticks(j) ticks(j)],[0 0.2],'k','linewidth',1)
 % end
 ticks = [ticks(1:2),ticks(end-3)];
 xticks(ticks)
 dd.XTickLabel = ["","","$\zeta$"];
 thisFig.Position = [351,156,916,643];
 dd.FontSize = 100;
 ylim([0 4])
 dd.TickLength = [0.01 0.1];
 dd.XAxis.TickLength(1) = 0.02;
 dd.XAxis.MajorTickChild.LineWidth = 3;
% dd.XAxis.MajorTickChild.
 set(gca, 'SortMethod', 'depth');
 set(gca ,'Layer', 'Top')
 
 % Adjust label position
 text(-0.07,-0.795735129068462,'$0$','Interpreter','latex','FontSize',100)
 text(0.04,-0.795735129068462,'$s_1$','Interpreter','latex','FontSize',100)
 orient(thisFig,'landscape')
 print(thisFig,'arcsineDistrPresent.pdf','-dpdf','-bestfit')
else
 orient(thisFig,'landscape')
 saveas(thisFig,'./../../paper_wwt/figures/arcsineDistr.pdf')
end




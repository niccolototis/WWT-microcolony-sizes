function [] = plotLastSlice(P,aSt) 
global spNames speciesPBE plotPVF
f = figure
x = aSt.AUT.pl.x2_all_TRUE_micronCube;
if plotPVF
 simul = aSt.AUT.pl.Marg_n_pvf_x2_RESC0to1; 
 if strcmp(P.initialPdf,'x2MarginalByKDEstimate') 
    ref_data = aSt.AUT.pl.n_pvf_data(1,:);
    ref_qH = aSt.AUT.pl.n_pvf_qH(1,:);
    ref_qL = aSt.AUT.pl.n_pvf_qL(1,:);
 end
else
 simul = aSt.AUT.pl.Marg_n_pdf_x2_RESC0to1;
 if strcmp(P.initialPdf,'x2MarginalByKDEstimate') 
    ref_data = aSt.AUT.pl.n_pdf_data(1,:);
    ref_qH = aSt.AUT.pl.n_pdf_qH(1,:);
    ref_qL = aSt.AUT.pl.n_pdf_qL(1,:);
 end
end

cmap = colormap(parula);
nCols = size(cmap,1);

hold on
fracCMap = 0.2;
thisRowCol = round(nCols*fracCMap);

fill([x(1) x x(end)],[0 ref_qH 0],cmap(thisRowCol,:),'FaceAlpha',0.6,'EdgeColor','none'); % Plots only the shapes with no stroke
fill([x(1) x x(end)],[0 ref_qL 0],cmap(thisRowCol,:),'FaceAlpha',1,'EdgeColor','none'); % Plots only the shapes with no stroke

% simulated
y = simul(end,:);
fracCMap = 0.9;
thisRowCol = round(nCols*fracCMap);
fill([x(1) x x(end)],[0 y 0],cmap(thisRowCol,:),'FaceAlpha',1,'EdgeColor','none');

plot(x,ref_data,'b','LineWidth', 1.5)

notshow = plot(x,ref_qL,'--b','LineWidth', 0.5);
set(get(get(notshow,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

notshow = plot(x,ref_qH,'--b','LineWidth', 0.5);
set(get(get(notshow,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

hold off


lgd = legend({'Experimental' 'Confidence intervals' 'Simulated'});
% lgd.FontSize = 15;
xlabel('Microcolony volume [\mum^3]','interpreter','Tex','fontsize',15);
ylabel('Density','fontsize',15) 
title('Steady-state PSD of autotrophs','FontSize',18) %select the name you want

% adjust
% - Ticks and grid
axes = gca;
axes.TickDir = 'out';
axes.YLim = [1400 1700];
axes.YTick = [1400 1500 1600 1700];         
axes.YAxis.Exponent = 3;

% - Axes labels position
f.OuterPosition(3) = f.OuterPosition(3)+100;
f.OuterPosition(4) = f.OuterPosition(4)+100;
set(axes,'OuterPosition',[0.1 0.1 0.9 0.9]);
axes.XLabel.Units = 'pixels';
axes.YLabel.Units = 'pixels';
axes.XLabel.Position(2) = axes.XLabel.Position(2)-2;
axes.YLabel.Position(1) = axes.YLabel.Position(1)-4;

% - Title
axes.Title.FontSize = 18;
axes.Title.Units = 'pixels';
axes.Title.Position(2) = axes.Title.Position(2)+5;


hold off
set(gcf,'Visible','on')


return

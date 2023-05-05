function [] = saveFigs(joyFigs,foldFig)
global IC_at labLambda improvStr cycleModStr plotOnlyFirstEndSlices plotForPresent plotOnlyEstimatedDistr
switch IC_at
 case 'SS_computed'
    lab = '_';
 case 'startupScenario'
    lab = '_startupScenario_';
end
if plotForPresent
 saveas(joyFigs.fig1,[foldFig 'joy' lab 'lambda_' labLambda improvStr cycleModStr '_presentaz.fig'])
 saveas(joyFigs.fig1,[foldFig 'joy' lab 'lambda_' labLambda improvStr cycleModStr '_presentaz.pdf'])
 saveas(joyFigs.fig1,[foldFig 'joy' lab 'lambda_' labLambda improvStr cycleModStr '_presentaz.jpeg'])
else
 saveas(joyFigs.fig1,[foldFig 'joy' lab 'lambda_' labLambda improvStr cycleModStr '.fig'])
 saveas(joyFigs.fig1,[foldFig 'joy' lab 'lambda_' labLambda improvStr cycleModStr '.pdf'])
end
if plotOnlyEstimatedDistr
 if plotForPresent
    saveas(joyFigs.fig2,[foldFig 'NDFest_presentaz' improvStr cycleModStr '.fig'])
    saveas(joyFigs.fig2,[foldFig 'NDFest_presentaz' improvStr cycleModStr '.pdf'])
    saveas(joyFigs.fig2,[foldFig 'NDFest_presentaz' improvStr cycleModStr '.jpeg'])
 else
    saveas(joyFigs.fig2,[foldFig 'NDFest' improvStr cycleModStr '.fig'])
    saveas(joyFigs.fig2,[foldFig 'NDFest' improvStr cycleModStr '.pdf'])
 end
end
if plotOnlyFirstEndSlices
 if plotForPresent
    saveas(joyFigs.fig3,[foldFig 'NDFend_' labLambda '_presentaz' improvStr cycleModStr '.fig'])
    saveas(joyFigs.fig3,[foldFig 'NDFend_' labLambda '_presentaz' improvStr cycleModStr '.pdf'])
    saveas(joyFigs.fig3,[foldFig 'NDFend_' labLambda '_presentaz' improvStr cycleModStr '.jpeg'])
 else
    saveas(joyFigs.fig3,[foldFig 'NDFend_' labLambda improvStr cycleModStr '.fig'])
    saveas(joyFigs.fig3,[foldFig 'NDFend_' labLambda improvStr cycleModStr '.pdf'])
 end
end
end

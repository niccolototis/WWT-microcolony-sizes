function [] = plotToSampleColors(start_iT)
% global noTrials
%     % use plot filled shapes with the right colors that I later sample
%     % in keynote and use the hexcolor codes to plot the lines in
%     % plot_ASM1seq_lowSRT. The problem is that there is no transparency for
%     % line objects in matlab    
%     colorbaseHET = [43 180 58]/255;
%     colorbaseAUT = hex2rgb('E5856E');
%     %
%     figure
%     hold on
%     for gT = start_iT:noTrials  
%         FaceAlpha = 0.4+gT*0.2;
%         aggt = fill([gT, gT, gT+1, gT+1],[0, 1, 1,0],colorbaseAUT,'FaceAlpha',FaceAlpha,'EdgeColor','none'); 
%         uistack(aggt,'top')
%         aggt = fill([gT, gT, gT+1, gT+1],[1, 2, 2, 1],colorbaseHET,'FaceAlpha',FaceAlpha,'EdgeColor','none'); 
%         uistack(aggt,'top')
% 
%     end
%     saveas(gcf,'./outputs/figures/ASM1seq/colors_lowSRT.pdf');
end

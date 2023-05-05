function [] = plotMovie(tpoints,x1_all_res,x2_all_res,Marg_x1_norm,Marg_x2_norm)
figure
for jj = 1:tpoints
 jj_s = strcat(['a',num2str(jj)]);
% plot(x1_all_res(jj,:), Marg_x1_norm.(jj_s), 'LineWidth', 2) 
% hold on
 plot(x2_all_res(jj,:), Marg_x2_norm.(jj_s), 'LineWidth', 2)  
% hold off
 legend('Marg x1')
 axis square
 title(jj_s)
 drawnow;
 pause(0.25);
end
return

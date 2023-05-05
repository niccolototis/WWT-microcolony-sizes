function [t_res_base,X_c_res_base,S_c_res_base,V_res_base] = runASM1seqbaselineParams(P,namePlot) 
global folderOutput T_sim
 % - temporarily change P.parSetup in this function
 parSetup_tmp = P.parSetup;
 P.parSetup = 'baseline';
 [t_res_base,X_c_res_base,S_c_res_base,V_res_base] = runASM1seqForOneParam(1,P);

 % - write to file   
 if P.writeASM1seqtoFile
    save([folderOutput,'resASM1seq_BASE_',num2str(T_sim),'.mat'],'t_res_base','X_c_res_base','S_c_res_base','V_res_base', '-v7.3');
 end
 
 if P.plotASM1seq
    figBase = figure;
    P = plotOneRunASM1seq(P,figBase,X_c_res_base(:,P.HET),'HET',namePlot,namePlot,1);
    P = plotOneRunASM1seq(P,figBase,X_c_res_base(:,P.AUT),'AUT',namePlot,namePlot,1);
 end
 P.parSetup = parSetup_tmp;   
return

function [P,figs_ASM1seq] = plot_ASM1seq_lowSRT(P,figs_ASM1seq)
global folderOutput T_sim cycleModStr improvStr iT 
% Initialize options, counters and flags
P = addToRecover(P);
P = settingsASM1seq(P); 

if P.plotASM1seq && ~P.plotPhases 
 if (isempty(figs_ASM1seq) && (iT ==-3 || iT==0))
    figs_ASM1seq.figHET = figure('visible','on','Position',[4,345,592,460]);
    figs_ASM1seq.figHET.Name = 'figHET';
    figs_ASM1seq.figAUT = figure('visible','on','Position',[421,345,592,460]);
    figs_ASM1seq.figAUT.Name = 'figAUT';
 end
end


 % - Run the ancillary scripts to have global variables defined for later usage
 P = settingsASM1seq(P); 
 P = defineKineticParams(P);
 
 % - I run setupSensAnSweep one time before initial
 % conditions because I need to have P.paramSetup defined. In
 % alternative put the definition of P.parSetup in generalSettings
 P = setupSensAnSweep(P);
 P = initializeSmatrix(P); 
 P = defineInletConcentrations(P);
 P = defineInitialConditions_t0(P);
 
 % - Load previously saved data
 % - Assign parameters for this run
 switch P.parSetup
   case 'chosenSetASM1seq'
        P = changeParam(P,1,P.parToChange);
        if ~isfile([folderOutput,'resASM1seq_',num2str(T_sim),improvStr,cycleModStr,'.mat'])
            error('Run the ASM1seq model first')
        else
            load([folderOutput,'resASM1seq_',num2str(T_sim),improvStr,cycleModStr,'.mat'],'t_res','X_c_res','S_c_res','V_res');
        end
        if P.plotSubs
            figure
            plot(S_c_res(:,P.sS))
            hold on
            plot(S_c_res(:,P.sNH4))
            plot(S_c_res(:,P.sND))
            plot(S_c_res(:,P.sNO))
            plot(S_c_res(:,P.sO2))
            plot(S_c_res(:,P.sI))
            legend({'sS','sNH4','sND','sNO','sO2','sI'})

            % hold on
            fprintf('SRT = %g, sS=%g sNH4=%g sND=%g sNO=%g sO2=%g sI=%g \n',P.SRT,mean(S_c_res(:,P.sS)),mean(S_c_res(:,P.sNH4)),mean(S_c_res(:,P.sND)),mean(S_c_res(:,P.sNO)),mean(S_c_res(:,P.sO2)),mean(S_c_res(:,P.sI)));
            cist = X_c_res(end,P.HET)/X_c_res(end,P.AUT);
            fprintf('SRT = %g, HA=%g \n',P.SRT,cist);
%             plot([P.tspallin(1),P.tspallin(end)],repmat(meanS,2,1))
            title(['SRT = ',num2str(P.SRT)])
            mu_ASM1seq_tall_HET = compute_mu(S_c_res,'HET',P);
            mu_ASM1seq_tall_AUT = compute_mu(S_c_res,'AUT',P); 
            fprintf('SRT = %g, mean mu_HET=%g, mean mu_AUT=%g \n',P.SRT,mean(mu_ASM1seq_tall_HET),mean(mu_ASM1seq_tall_AUT));
            figure
            plot(mu_ASM1seq_tall_HET)
            hold on
            plot(mu_ASM1seq_tall_AUT)
            title(['mu for SRT = ',num2str(P.SRT)])
        end
        if P.plotASM1seq
            if P.plotPhases
                error('non prevedo di plottare questo')
            else
                figs_ASM1seq = plotOneRunASM1seq_lowSRT(figs_ASM1seq,P,X_c_res(:,P.HET),'HET');
                figs_ASM1seq = plotOneRunASM1seq_lowSRT(figs_ASM1seq,P,X_c_res(:,P.AUT),'AUT');
                figs_ASM1seq = adjustPlotAndSave_lowSRT(figs_ASM1seq,P,'HET');
                figs_ASM1seq = adjustPlotAndSave_lowSRT(figs_ASM1seq,P,'AUT');
            end
        end
    otherwise
        error('should be P.parSetup set to chosenSetASM1seq')
 end

return

function P = addToRecover(P)
% Here I assign the fields of P to make this file able to call the
% dependent funtions present in run..._bis
P.modelStructure = 'actualSystem'; 
P.modelICPars = 'actualSystem';
P.modelKinPars = 'actualSystem';%'ASM1'
P.ICisSS = true; 
P.ICisSCells = false;
return









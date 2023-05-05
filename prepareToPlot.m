function aSt = prepareToPlot(aSt,aSt_options,varargin)
global speciesPBE spNames T_sim comparison
% compute mean and standard deviation using expected values
warning ('off','all');
if nargin>2
 compPD = varargin{1};
end
for aSpPBE = speciesPBE
 sp = char(spNames{aSpPBE});
 Ss = aSt.(sp);
 x2_ind_toPlot = aSt_options.(sp).x2_ind_toPlot;
 x2_all = Ss.x2_all_TRUE_micronCube(:)';
 x2_b = aSt_options.(sp).x2_b_TRUE_micronCube(:)';
 x2_lb = x2_b(1:end-1);
 x2_ub = x2_b(2:end);
 % x2_widths = diff(x2_b);
 x2min = min(x2_b);
 x2max = max(x2_b);
 nTpoints = size(Ss.Marg_N_x2_RESC0to1,1);  
 % NB this is to improve the estimate
 xall = sort([x2_all x2_b]);   
 
 %% data
 % Take compPD -> mean_model = mean_data for t=0 when IC_at=SS_computed
 [mean_data,q1_data,q3_data] = deal(ones(nTpoints,1));
 folderOutput_tmp = ['outputs/data/Tsim_',num2str(T_sim),'_SS/'];
 fold_S_compPD_tmp = [folderOutput_tmp 'S_compPD/'];
 load([fold_S_compPD_tmp 'S_compPD_t0.mat'],'compPD_t0')
 compPD_t0 = truncate(compPD_t0,x2min,x2max);
 quantiles = cdf(compPD_t0,xall);
 [~,indmean] = min(abs(quantiles-0.5));
 mean_model_t0 = xall(indmean);
 mean_data = mean_data.*mean_model_t0;
 [~,ind1] = min(abs(quantiles-0.15));
 q1_model_t0 = xall(ind1);
 [~,ind3] = min(abs(quantiles-0.75));
 q3_model_t0 = xall(ind3);
 q1_data = q1_data.*q1_model_t0;
 q3_data = q3_data.*q3_model_t0;
 
 %% Model
 [mean_model,q1_model,q3_model] = deal(zeros(nTpoints,1));
 % Truncate the probability distributions!
 for oneTpoint = 1:nTpoints
    aTp = ['t' num2str(oneTpoint)];
    compPD.(aTp) = truncate(compPD.(aTp),x2min,x2max);
    quantiles = cdf(compPD.(aTp),xall);
    [~,indmean] = min(abs(quantiles-0.5));
    mean_model(oneTpoint) = xall(indmean);
    [~,ind1] = min(abs(quantiles-0.15));
    q1_model(oneTpoint) = xall(ind1);
    [~,ind3] = min(abs(quantiles-0.75));
    q3_model(oneTpoint) = xall(ind3);
 end  
 
 % for the exp data I concatenate for each time point always the same
 % vector
 N_data_tall = repmat(Ss.N_data(:)',nTpoints,1);
 N_qL_tall = repmat(Ss.N_qL(:)',nTpoints,1);
 N_qH_tall = repmat(Ss.N_qH(:)',nTpoints,1);
 
 
 %% - Put in struct
 pl.nTpoints = nTpoints;
 pl.x2_all_TRUE_micronCube = x2_all(x2_ind_toPlot);
 pl.x2_b_TRUE_micronCube = x2_b([x2_ind_toPlot x2_ind_toPlot(end)+1]);
 % - N
 pl.Marg_N_x2_RESC0to1 = Ss.Marg_N_x2_RESC0to1(:,x2_ind_toPlot);
 pl.N_data = N_data_tall(:,x2_ind_toPlot);
 pl.N_qL = N_qL_tall(:,x2_ind_toPlot);
 pl.N_qH = N_qH_tall(:,x2_ind_toPlot);
 % - pdf
 pl.Marg_n_pdf_x2_RESC0to1 = Ss.Marg_n_pdf_x2_RESC0to1(:,x2_ind_toPlot);
 pl.n_pdf_data = repmat(Ss.n_pdf_data(x2_ind_toPlot)',nTpoints,1);
 pl.n_pdf_qL = repmat(Ss.n_pdf_qL(x2_ind_toPlot)',nTpoints,1);
 pl.n_pdf_qH = repmat(Ss.n_pdf_qH(x2_ind_toPlot)',nTpoints,1);
 % - pvf
 pl.Marg_n_pvf_x2_RESC0to1 = Ss.Marg_n_pvf_x2_RESC0to1(:,x2_ind_toPlot);
 pl.n_pvf_data = repmat(Ss.n_pvf_data(x2_ind_toPlot)',nTpoints,1);
 pl.n_pvf_qL = repmat(Ss.n_pvf_qL(x2_ind_toPlot)',nTpoints,1);
 pl.n_pvf_qH = repmat(Ss.n_pvf_qH(x2_ind_toPlot)',nTpoints,1);
 % - mean
 pl.mean_model = mean_model;
 pl.mean_data = mean_data;
 % - 1 annd 3 quartiles
 pl.q1_model = q1_model;
 pl.q3_model = q3_model;
 pl.q1_data = q1_data;
 pl.q3_data = q3_data;

 % -
 pl.t_all_PBE = aSt.(sp).t_all_PBE;   
 aSt.(sp).pl = pl;
end 
warning ('on','all'); 
end

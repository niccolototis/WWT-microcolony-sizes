function [joyFigs,P] = prepareJoyPlot_lowerSRT(P,aSt,varargin)
global spNames speciesPBE plotWhat IC_at iT

for aSpPBE = speciesPBE
 sp = char(spNames{aSpPBE});
 pl = aSt.(sp).pl;
 if nargin>2
    aSt_SS = varargin{1};
    if ~isempty(aSt_SS)
        pl_SS = aSt_SS.(sp).pl;
    end
 end
 if nargin>3
    joyFigs = varargin{2};
 end

 switch IC_at
    case 'SS_computed'
        %%
        % - Select 11 frames
        nSlices = 11;
        idSlices = round(linspace(1,P.ntpts_plot_PBE,nSlices));
        yLabels = aSt.AUT.t_all_PBE(idSlices);

        % s1 is model (impo to transpose slices are given in columns)
        x_struct.s1 = pl.x2_all_TRUE_micronCube';
        switch plotWhat
            case 'pvf'
                data_struct.s1.model = pl.Marg_n_pvf_x2_RESC0to1(idSlices,:)'; 
            case 'N'
                data_struct.s1.model = pl.Marg_N_x2_RESC0to1(idSlices,:)'; 
            case 'pdf'
                data_struct.s1.model = pl.Marg_n_pdf_x2_RESC0to1(idSlices,:)'; 
        end
        data_struct.s1.mean_model = pl.mean_model(idSlices)';
        data_struct.s1.q1_model = pl.q1_model(idSlices)';
        data_struct.s1.q3_model = pl.q3_model(idSlices)';
        showSlices_struct.s1 = 2:nSlices; % exclude the first!
        plot_type_struct.s1 = 'filled';

        Facecolor_struct.s1 = linspace(0.2,0.9,nSlices);
        StrokeColor_struct.s1 = [];

        if strcmp(P.initialNDF,'x2MarginalByKDEstimate')   
            
            % s2 is data
            x_struct.s2 = pl.x2_all_TRUE_micronCube';
            switch plotWhat
                case 'pvf'
                    data_struct.s2.data = pl.n_pvf_data(idSlices,:)'; 
                case 'N'
                    data_struct.s2.data = pl.N_data(idSlices,:)'; 
                case 'pdf'
                    data_struct.s2.data = pl.n_pdf_data(idSlices,:)'; 
            end
            data_struct.s2.mean_data = pl.mean_data(idSlices)';
            data_struct.s2.q1_data = pl.q1_data(idSlices)';
            data_struct.s2.q3_data = pl.q3_data(idSlices)';
            if iT ==0
                showSlices_struct.s2 = [1];                
            else
                showSlices_struct.s2 = [];                
            end
            plot_type_struct.s2 = 'filled';
            Facecolor_struct.s2 = 0.2*ones(1,nSlices);
            StrokeColor_struct.s2 = 'b';

            % s3 is qL
            x_struct.s3 = pl.x2_all_TRUE_micronCube';
            switch plotWhat
                case 'pvf'
                    data_struct.s3.qL = pl.n_pvf_qL(idSlices,:)'; 
                case 'N'
                    data_struct.s3.qL = pl.N_qL(idSlices,:)';
                case 'pdf'
                    data_struct.s3.qL = pl.n_pdf_qL(idSlices,:)'; 
            end
            if iT ==0
                showSlices_struct.s3 = [1];                
            else
                showSlices_struct.s3 = [];                
            end
            plot_type_struct.s3 = 'filled';
            Facecolor_struct.s3 = 0.2*ones(1,nSlices);
            StrokeColor_struct.s3 = [];

            % s4 is qH
            x_struct.s4 = pl.x2_all_TRUE_micronCube';
            switch plotWhat
                case 'pvf'
                    data_struct.s4.qH = pl.n_pvf_qH(idSlices,:)'; 
                case 'N'
                    data_struct.s4.qH = pl.N_qH(idSlices,:)'; 
                case 'pdf'
                    data_struct.s4.qH = pl.n_pdf_qH(idSlices,:)'; 
            end
            if iT ==0
                showSlices_struct.s4 = [1];                
            else
                showSlices_struct.s4 = [];                
            end
            plot_type_struct.s4 = 'filled';
            Facecolor_struct.s4 = 0.2*ones(1,nSlices);
            StrokeColor_struct.s4 = [];
        end
        overlapMethod = 'constant';         
        
    case 'startupScenario' 
        %%
        % - Select 11 frames
        idSlices = [1 2 3 4 7 11 15 21 31 41 51 61 71 81 91 101];
        nSlices = length(idSlices);
        yLabels = aSt.AUT.t_all_PBE(idSlices);

        % s1 is model (impo to transpose slices are given in columns)
        x_struct.s1 = pl.x2_all_TRUE_micronCube';
        switch plotWhat
            case 'pvf'
                data_struct.s1.model = pl.Marg_n_pvf_x2_RESC0to1(idSlices,:)'; 
            case 'N'
                data_struct.s1.model = pl.Marg_N_x2_RESC0to1(idSlices,:)'; 
            case 'pdf'
                data_struct.s1.model = pl.Marg_n_pdf_x2_RESC0to1(idSlices,:)'; 
        end
        data_struct.s1.mean_model = pl.mean_model(idSlices)';
        data_struct.s1.q1_model = pl.q1_model(idSlices)';
        data_struct.s1.q3_model = pl.q3_model(idSlices)';
        showSlices_struct.s1 = 1:nSlices; % exclude the first!
        plot_type_struct.s1 = 'filled';
        Facecolor_struct.s1 = linspace(0.2,0.9,nSlices);
        StrokeColor_struct.s1 = [];
        
        % s5 is  SS_computed model for CFR with low_biomass
        x_struct.s5 = pl_SS.x2_all_TRUE_micronCube';
        switch plotWhat
            case 'pvf'
                data_struct.s5.model = pl_SS.Marg_n_pvf_x2_RESC0to1(idSlices,:)'; 
            case 'N'
                data_struct.s5.model = pl_SS.Marg_N_x2_RESC0to1(idSlices,:)'; 
            case 'pdf'
                data_struct.s5.model = pl_SS.Marg_n_pdf_x2_RESC0to1(idSlices,:)'; 
        end
        data_struct.s5.mean_model = pl_SS.mean_model(idSlices)';
        data_struct.s5.q1_model = pl_SS.q1_model(idSlices)';
        data_struct.s5.q3_model = pl_SS.q3_model(idSlices)';
        showSlices_struct.s5 = 1:nSlices; % exclude the first!
        plot_type_struct.s5 = 'filled';
        Facecolor_struct.s5 = linspace(0.2,0.9,nSlices);
        StrokeColor_struct.s5 = [];

        if strcmp(P.initialNDF,'x2MarginalByKDEstimate')  
            
            % s2 is data
            x_struct.s2 = pl.x2_all_TRUE_micronCube';
            switch plotWhat
                case 'pvf'
                    data_struct.s2.data = pl.n_pvf_data(idSlices,:)'; 
                case 'N'
                    data_struct.s2.data = pl.N_data(idSlices,:)'; 
                case 'pdf'
                    data_struct.s2.data = pl.n_pdf_data(idSlices,:)'; 
            end
            data_struct.s2.mean_data = pl.mean_data(idSlices)';
            data_struct.s2.q1_data = pl.q1_data(idSlices)';
            data_struct.s2.q3_data = pl.q3_data(idSlices)';
            showSlices_struct.s2 = [];
            plot_type_struct.s2 = 'filled';
            Facecolor_struct.s2 = 0.2*ones(1,nSlices);
            StrokeColor_struct.s2 = 'b';

            % s3 is qL
            x_struct.s3 = pl.x2_all_TRUE_micronCube';
            switch plotWhat
                case 'pvf'
                    data_struct.s3.qL = pl.n_pvf_qL(idSlices,:)'; 
                case 'N'
                    data_struct.s3.qL = pl.N_qL(idSlices,:)';
                case 'pdf'
                    data_struct.s3.qL = pl.n_pdf_qL(idSlices,:)'; 
            end
            showSlices_struct.s3 = [];
            plot_type_struct.s3 = 'filled';
            Facecolor_struct.s3 = 0.2*ones(1,nSlices);
            StrokeColor_struct.s3 = [];

            % s4 is qH
            x_struct.s4 = pl.x2_all_TRUE_micronCube';
            switch plotWhat
                case 'pvf'
                    data_struct.s4.qH = pl.n_pvf_qH(idSlices,:)'; 
                case 'N'
                    data_struct.s4.qH = pl.N_qH(idSlices,:)'; 
                case 'pdf'
                    data_struct.s4.qH = pl.n_pdf_qH(idSlices,:)'; 
            end   
            showSlices_struct.s4 = [];
            plot_type_struct.s4 = 'filled';
            Facecolor_struct.s4 = 0.2*ones(1,nSlices);
            StrokeColor_struct.s4 = [];
        end
        overlapMethod = 'variable';
 end
 x_struct = orderfields(x_struct);
 data_struct = orderfields(data_struct);
 plot_type_struct = orderfields(plot_type_struct);
 showSlices_struct = orderfields(showSlices_struct);
 Facecolor_struct = orderfields(Facecolor_struct);
 StrokeColor_struct = orderfields(StrokeColor_struct);

 [joyFigs,P] = joyPlot_lowerSRT(P,x_struct,data_struct,plot_type_struct,showSlices_struct,yLabels,-0.3,...
    'variable',false,'FaceColor',Facecolor_struct,'StrokeColor',StrokeColor_struct,'LineColor','black',...
    'overlapMethod',overlapMethod,'VLines',true,'joyFigs',joyFigs);    
end
return


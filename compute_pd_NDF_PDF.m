function [varargout] = compute_pd_NDF_PDF(sampledSizes,x2_b,xforPDF,width,outWhat)
 % d] Estimating NDF at grid points x2_all. Note that x2_all is not an
 % input argument, only x2_b is required 
 %outWhat = {'pd','NDF','PDF'}
 if ~(iscell(outWhat) && ischar(outWhat{1})) 
    error('wrong input, outWhat needs to be a cell')
 end
 sampledSizes = sort(sampledSizes);
 pd = fitdist(sampledSizes,'Kernel','Kernel','normal','Width',width);
 for i = 1:length(outWhat)
    oneOut = outWhat{i};
    switch oneOut
        case 'pd'
            varargout{i} = pd;
        case 'NDF'
            x2_lb = x2_b(1:end-1);
            x2_ub = x2_b(2:end);
            expCDF_at_x2_lb = cdf(pd,x2_lb);
            expCDF_at_x2_ub = cdf(pd,x2_ub);
            N = expCDF_at_x2_ub-expCDF_at_x2_lb;
            varargout{i} = N;
        case 'PDF'
            PDF = pdf(pd,xforPDF);
            varargout{i} = PDF;
    end
 end
return

function eCDF = compute_eCDF(data_sizes,whereToEvalCDF)
 % x_all are the 
 n = length(data_sizes);
 sizes_M = repmat(data_sizes(:)',length(whereToEvalCDF),1);
 bool_M = (sizes_M<=whereToEvalCDF(:));
 eCDF = 1/n*sum(double(bool_M),2);
end

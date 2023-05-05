function [nTpoints_plot,varargout] = shrinkSize(nTpoints_all,varargin)
if nTpoints_all>1000
 Q = floor(nTpoints_all/100);
 idxs = 1:Q:nTpoints_all;
 idxs = [idxs nTpoints_all];
 nTpoints_plot = length(idxs);
 for i = 1:(nargin-1)
    if isvector(varargin{i})
        varargout{i} = varargin{i}(idxs);
    else
        varargout{i} = varargin{i}(idxs,:);
    end
 end
else 
 nTpoints_plot = nTpoints_all;
 varargout = varargin;
end
return

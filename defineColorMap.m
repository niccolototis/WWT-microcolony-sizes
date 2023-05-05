function P = defineColorMap(P)
%-- Set color scheme --%
colormap hsv

P.allColors.baseLiterature = [65,65,65]/255;

% color map with green [0 1 0] and blue [0 0 1] on the edges
C.HET = ['00D637'; 'E1FB70']; % green to boh
C.HET = hex2rgb(C.HET);
C.AUT = ['CFAAA8';'E1D480']; % red to boh
C.AUT = hex2rgb(C.AUT);

P.allColors.baseHET = [43 180 58]/255;
P.allColors.baseAUT = hex2rgb('E5856E');

P.allColors.phaseColors = ...
{[179 237 255]/255,...
[203 69 160]/255,...
[146 146 146]/255,...
[229 133 110]/255,...
[100 168 255]/255};

for cc = {'HET','AUT'}
 P.allColors.(cc{1}) = nan(P.nSweeps,3); 
 P.allColors.(cc{1})(:,1) = linspace(C.(cc{1})(1,1),C.(cc{1})(2,1),P.nSweeps)';
 P.allColors.(cc{1})(:,2) = linspace(C.(cc{1})(1,2),C.(cc{1})(2,2),P.nSweeps)';
 P.allColors.(cc{1})(:,3) = linspace(C.(cc{1})(1,3),C.(cc{1})(2,3),P.nSweeps)';
end
end

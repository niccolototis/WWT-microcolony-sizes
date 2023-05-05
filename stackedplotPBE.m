function [] = stackedplotPBE(x_s,data_s,plot_type_s,showSlices_s,offset,varargin)
%% Parse the inputs
% Initialize the input parser
p = inputParser;
% Set defaults
dOverlapMethod = 'constant';
dReverse = false;
dFaceColor = 'w';
dFaceAlpha = 1;
dLineColor = 'k';
dLineWidth = 1.5;
dStrokeColor = 'w';
dStrokeWidth = 3;
dVLines = 0;
dVLinesColor = '--k';
dVLinesWidth = 1;
didxFramesToAdd = [];
% Validation function
expectedMethods = {'variable','constant'};
validOverlapMethod = @(x) isnumeric(x) ||...
 (ischar(x) && any(validatestring(x,expectedMethods)));
% Add required
addRequired(p,'data_s',@isstruct);
addRequired(p,'x_s',@isstruct);
addRequired(p,'plot_type_s',@isstruct);
addRequired(p,'showSlices_s',@isstruct);
addRequired(p,'offset',@isnumeric);

% Add optional parameters
addOptional(p,'overlapMethod',dOverlapMethod,validOverlapMethod)
addOptional(p,'reverse',dReverse,@islogical)
% Add pair parameters
addParameter(p,'FaceColor',dFaceColor)
addParameter(p,'StrokeColor',dStrokeColor)
addParameter(p,'FaceAlpha',dFaceAlpha,@isnumeric)
addParameter(p,'LineColor',dLineColor)
addParameter(p,'LineWidth',dLineWidth,@isnumeric)
addParameter(p,'StrokeWidth',dStrokeWidth,@isnumeric)
addParameter(p,'VLines',dVLines,@isnumeric)
addParameter(p,'VLinesColor',dVLinesColor)
addParameter(p,'VLinesWidth',dVLinesWidth,@isnumeric)
addParameter(p,'idxFramesToAdd',didxFramesToAdd,@isnumeric)
% Parse
parse(p,x_s,data_s,plot_type_s,showSlices_s,offset,varargin{:})
overlapMethod = p.Results.overlapMethod;
reverse = p.Results.reverse;
FaceColor_s = p.Results.FaceColor;
StrokeColor_s = p.Results.StrokeColor;
FaceAlpha = p.Results.FaceAlpha;
LineColor = p.Results.LineColor;
LineWidth = p.Results.LineWidth;
StrokeWidth = p.Results.StrokeWidth;
VLines = p.Results.VLines;
VLinesColor = p.Results.VLinesColor;
VLinesWidth = p.Results.VLinesWidth;
idxFramesToAdd = p.Results.idxFramesToAdd;

%% Data preparation
% If the user wants to return the handles, so be it
nout = max(nargout,1) - 1;
cnt = 0;
fn_x = fieldnames(x_s);
fn_data = fieldnames(data_s);
fn_plot_type = fieldnames(plot_type_s);
fn_SS = fieldnames(showSlices_s);
if ~isequal(fn_x,fn_data,fn_plot_type,fn_SS)
 ME = MException('data_struc and x_struc should contain the same fields');
 throw(ME);
end
if ~isstruct(FaceColor_s)
 FaceColor_in = FaceColor_s;
 FaceColor_s = [];
 for iss = 1:length(fn_x)
    ss = char(fn_x(iss));
    FaceColor_s.(ss) = FaceColor_in;
 end
end
if ~isstruct(StrokeColor_s)
 StrokeColor_in = StrokeColor_s;
 StrokeColor_s = [];
 for iss = 1:length(fn_x)
    ss = char(fn_x(iss));
    StrokeColor_s.(ss) = StrokeColor_in;
 end
end
 
%% Start plotting from the dataset most in the back to the front one
for nS = numel(fn_data):-1:1 
 ss = fn_data{nS};
 x = x_s.(ss);
 plot_type = plot_type_s.(ss);
 if ~isnumeric(data_s.(ss)) || ~isnumeric(x)
    ME = MException('Numerical values are needed');
    throw(ME);
 end    
end

s = stackedplot(outdoors,'MarkerFaceColor','r');
s.LineWidth = 2;
s.LineProperties
s.LineProperties(2).PlotType = 'scatter';
s.LineProperties(3).PlotType = 'stairs';
s.AxesProperties
end

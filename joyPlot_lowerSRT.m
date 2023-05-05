function [joyFigs,P] = joyPlot_lowerSRT(P,x_s,data_s,plot_type_s,showSlices_s,yTickLabels,offset,varargin)
global  IC_at noTrials labLambda improvStr cycleModStr  figsVisible iT
%JOYPLOT Plot your data in a ridgeline representation
% JOYPLOT(DATA,X,OFFSET) The array DATA should have a size m by n, where 
% n is the number ofdatasets to plot and m is the sampling number. X is 
% a vector containing the x-coordinates. OFFSET is a scalar that 
% determines the displacement between plots. 
%
% JOYPLOT(DATA,X,OFFSET,OVERLAPMETHOD) If the optional parameter
% OVERLAPMETHOD is set to 'variable', OFFSET is a the overlap between 
% datasets in percent. By default, OVERLAPMETHOD is set to 'constant'.
%
% JOYPLOT(DATA,X,OFFSET,OVERLAPMETHOD,REVERSE) If REVERSE is set to
% 'true', the first row of DATA will be plotted at the bottom and the
% last at the bottom. 
%
% JOYPLOT(___,Name,Value) Specifies patch and line properties using one
% or more 'Name,Value' pairs. For a list see Properties below.
%
% [hf,hl] = JOYPLOT(___) Returns the patch and line handles.
%
% [hs,hf,hl] = JOYPLOT(___) If the stroke is requested, its handle can be
% returned as well.
%
% [hs,hf,hl,hvl] = JOYPLOT(___) If additional vertical lines are plotted,
% their handles can be returned too.
%
% The function accepts the following Name,Value pairs:
%
% - FaceColor: works just like the color in the PATCH function. It
% additionally accepts the input 'position', which colors the faces
% according to the position in the y-axis.
% - FaceAlpha: works just like 'FaceAlpha' in the function PATCH.
% - LineColor: works just like the color of function PLOT, with the
% addition that it can be set to 'none' to delete the line.
% - LineWidth: works just like 'LineWidth' in the function PLOT.
% - StrokeColor: adds a stroke to each dataset with the specified color.
% Here, the color works just like in the funciton PATCH.
% - StrokeWidth: adds a stroke to each dataset with the specified width.
% Again, this works like in the function PATCH.
% - VLines: adds vertical lines in each dataset corresponding to the
% x-coordinates provided. This is useful to plot median, mean or mode
% values for example.
% - VLinesColor: works just like the color of function PLOT.
% - VLinesWidth: works just like 'LineWidth' in the function PLOT.
%
% EXAMPLE
% The classic joy plot from the cover of Joy Division's Unknown Pleasures 
% by Peter Saville.
% % Save the data from the web
% filename = websave('pulsar.csv',['https://gist.githubusercontent.com/borgar/',...
% '31c1e476b8e92a11d7e9/raw/0fae97dab6830ecee185a63c1cee0008f6778ff6/',...
% 'pulsar.csv']);
% % Import it into MATLAB(R)
% data = readmatrix(filename);
% x = linspace(0,93,size(data,2));
% % Create the figure and show the magic
% figure
% joyPlot(data',x,4)
% set(gcf,'position',[500,100,560,680])
% set(gca,'Visible','off', 'box','off','XTick',[],'YTick',[])
%
% See also PLOT, PATCH, FILL.
%
% Author: Santiago Benito, Ruhr-Universität Bochum
% Contact: santiago.benito@rub.de
% Last modified: 23.04.2020
% Version: v0.1
% ----------------------------------------------------------------- @NT the
% datasets in input need to be given in the right order of how the plots
% need to be superimposed from the most in front to the most in the back


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
dVLinesColor = ':k';
dVLinesWidth = 1;
didxFramesToAdd = [];
djoyFigs = [];
% Validation function
expectedMethods = {'variable','constant'};
validOverlapMethod = @(x) isnumeric(x) ||...
 (ischar(x) && any(validatestring(x,expectedMethods)));
% Add required
addRequired(p,'data_s',@isstruct);
addRequired(p,'x_s',@isstruct);
addRequired(p,'plot_type_s',@isstruct);
addRequired(p,'showSlices_s',@isstruct);
addRequired(p,'yLabels',@isnumeric);
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
addParameter(p,'VLines',dVLines,@boolean)
addParameter(p,'VLinesColor',dVLinesColor)
addParameter(p,'VLinesWidth',dVLinesWidth,@isnumeric)
addParameter(p,'idxFramesToAdd',didxFramesToAdd,@isnumeric)
addParameter(p,'joyFigs',djoyFigs)
% Parse
parse(p,x_s,data_s,plot_type_s,showSlices_s,yTickLabels,offset,varargin{:})
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
joyFigs = p.Results.joyFigs;

%%
if isempty(joyFigs)
 fig1 = figure('Visible',figsVisible,'Position', [0, 100, 500, 800]);
else
 if isfield(joyFigs,'fig1')
    fig1 = figure(joyFigs.fig1);
 else
    error('this field should already present from having plotted the baseline results with parametersEffectOnMSD = false. That needs the first figure to be plotted, then the others with parametersEffectOnMSD=true can be added on top')
 end
end
hold on

%% Data preparation
% If the user wants to return the handles, so be it
nout = max(nargout,1) - 1;
cnt = 0;

% @NT inserted for loop here
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
%% Get the y-positions from the first dataset (the most in front)
data_tmp = data_s.s1.model;
[m_tmp,n_tmp] = size(data_tmp);
if ~reverse 
 data_tmp = fliplr(data_tmp);
end

% Defining the space between 
maxHeightExpData = max(max(data_s.s2.data,[],1));
ypos = cumsum(ones(1,size(data_tmp,2))*maxHeightExpData)*(1-offset);
ypos = [0, ypos]; ypos(end) = [];
switch overlapMethod
 case 'variable'
    %     ypos = cumsum(max(data_tmp,[],1))*(1-offset);
    %     ypos = [0, ypos]; ypos(end) = [];
    if min(ypos)~= 0
        error('should be zero')
    end
    maxYpos = max(ypos);
    % rescale
    ypos = yTickLabels(:)';
    ypos = (ypos/max(ypos))*maxYpos;
 case'constant'

    %     ypos = 0:offset:offset*(n_tmp-1); 
end

%ypos sets the y-coord that tells me the position on the vertical axis
%where I draw the orizontal lines of all the x-axes one for each frame
locationHorizLines = repmat(ypos,m_tmp,1);
locationHorizLines = max(max(locationHorizLines))-locationHorizLines;
locationHorizLines = fliplr(locationHorizLines);

% if strcmp(overlapMethod,'variable')
% howManyTimesDistribHigherThanSpacingInbetweenSlices = 5;
% % Rescale locationHorizLines to avoid that the curves representing the
% % data become too small wrt the distance among them
% minSepHeigh = min(diff(unique(locationHorizLines)));
% % I take one (cut) slice of data
% oneData = data_s.s2.data(:,1);
% % find the max (top of the distribution)
% maxHeightData = max(oneData);
% % Adjust locationHorizLines so that maxHeightData is twice as much as
% % the min distance among slices
% minSepHeigh_desired = maxHeightData/5;
% multFactor = minSepHeigh_desired/minSepHeigh;
% locationHorizLines = locationHorizLines*multFactor;
% end

%% Start plotting from the dataset most in the back to the front one
for nS = numel(fn_data):-1:1 
 ss = fn_data{nS};
 x = x_s.(ss);
 plot_type = plot_type_s.(ss);
 
 % Format x as a column vector
 x = x(:);
 
 % Extract some important parameters
 mini = min(x);
 maxi = max(x);
 
 switch ss
    case 's1'
        thisData = data_s.s1.model;
    case 's2'
        thisData = data_s.s2.data;
    case 's3'
        thisData = data_s.s3.qL;
    case 's4'
        thisData = data_s.s4.qH;
    case 's5'
        thisData = data_s.s5.model;
 end
    
 
 [m,n] = size(thisData);
 if m~= m_tmp || n~=n_tmp
    ME = MException('the dataset for the model, for the real data and for the confidence intervals need to have the same dimensions');
    throw(ME);
 end
 % Apply the reverse option if needed
 if ~reverse
    thisData = fliplr(thisData);
 end

 % Prepare the data for patch plotting
 %x is a column. I double the min and max values of x at the extremes
 %and I create n column equal to x, each column is the x-coord of one
 %frame
 %isequal(sort(X(2:end-1,1)),sort(x))
 
 Y.(ss) = thisData + locationHorizLines;
 switch plot_type
    case 'filled'
        X.(ss) = repmat([mini;x;maxi],1,n);
 %        Y.(ss) = fliplr([min(Y.(ss),[],1);Y.(ss);min(Y.(ss),[],1)]);
        Y.(ss) = fliplr([locationHorizLines(1,:);Y.(ss);locationHorizLines(1,:)]);
    case 'line'
        X.(ss) = repmat(x,1,n);
        Y.(ss) = fliplr(Y.(ss));
 end

 % Check if a colormap is to be imployed
 if strcmp(FaceColor_s.(ss),'position')
    FaceColor_s.(ss) = ypos;   
 end
 % Check if the stroke needs to be used
 stroke.(ss) = true;
% if sum(strcmp(varargin,'StrokeColor') + ...
%         strcmp(varargin,'StrokeWidth')) > 0
%     stroke.(ss) = true;
% end
 if isempty(StrokeColor_s.(ss)) || strcmp(StrokeColor_s.(ss),'w')
    stroke.(ss) = false;
 end
 %% Patch & Plot
 % Get the current axes, if given
 ax = gca;
 % Store the number of existing children
 numelAx = numel(ax.Children);
end

cmap = colormap(parula);
nCols = size(cmap,1);
hold on

%% - Plot mean of the experiment below everything
if iT ==0
 % yv = interp1(x,data_s.s2.data(:,aSlice),data_s.s2.mean_data(aSlice)); 
 XV = [data_s.s2.mean_data(1) data_s.s2.mean_data(1)];
 YV = [0 max(locationHorizLines(1,:))];
 coldata = 'r';
 if isempty(XV(XV<mini|XV>maxi))
    plot(XV,YV,'-','color',coldata,'LineWidth',VLinesWidth);
 end
 XV = [data_s.s2.q1_data(1) data_s.s2.q1_data(1)];
 YV = [0 max(locationHorizLines(1,:))];
 if isempty(XV(XV<mini|XV>maxi))
    plot(XV,YV,':','color',coldata,'LineWidth',VLinesWidth*5/3);
 end
 XV = [data_s.s2.q3_data(1) data_s.s2.q3_data(1)];
 YV = [0 max(locationHorizLines(1,:))];
 if isempty(XV(XV<mini|XV>maxi))
    plot(XV,YV,':','color',coldata,'LineWidth',VLinesWidth*5/3);
 end
end
%% - Plot model
for aSlice = 1:n
 for nS = 4:-1:1 
    ss = fn_data{nS};
    plot_type = plot_type_s.(ss);
    StrokeColor = StrokeColor_s.(ss);
    fracCMap = FaceColor_s.(ss)(aSlice);
    thisRowCol = round(nCols*fracCMap);
    % Plot the rest of the stuff
    if ismember(aSlice,showSlices_s.(ss)) 
        if strcmp(ss,'s2') && iT ==0  
            if~strcmp(LineColor,'none')
                if strcmp(LineColor,'darker')
                    tmpfracCMap = fracCMap-0.15;
                    thisRowCol = round(nCols*tmpfracCMap);
                    col = cmap(thisRowCol,:);
                    plot(X.(ss)(2:end-1,aSlice),Y.(ss)(2:end-1,aSlice),'Color',col,'LineWidth',LineWidth);
                else
                    if aSlice ==1
                        col = 'k';
                    else
                        col = StrokeColor_s.(ss);
                    end
                    plot(X.(ss)(:,aSlice),Y.(ss)(:,aSlice),col,'Linewidth',0.5); % Plots only the shapes with no stroke
                end
            end
        else
            switch plot_type
                case 'filled'
                    storeFaceAlpha = FaceAlpha;
                    switch ss
                        case 's1'
                            col = cmap(thisRowCol,:);
                            FaceAlpha = 1;
                            switch iT
                                case 0
                                    %FaceAlpha = 0.4;
                                    %col = cmap(thisRowCol,:);
                                    col = hex2rgb('32D644');
                                case 1
                                    %FaceAlpha = 0.6;
                                    %col = hex2rgb('FD6039');
                                    col = hex2rgb('C06A72');
                                case 2
                                    %FaceAlpha = 0.8;
                                    col = hex2rgb('7888C0');
                                case 3
                                    %FaceAlpha = 1;
                                    col = hex2rgb('C09954');
                            end
                        case 's2'
                            % already plotted
                        case 's3'
                            col = [0 0 0];
                            FaceAlpha = 0.5;
                        case 's4'
                            col = [0 0 0];
                            FaceAlpha = 0.2;
                    end
                    pp.(ss) = fill(X.(ss)(:,aSlice),Y.(ss)(:,aSlice),col,'FaceAlpha',FaceAlpha,'EdgeColor','none'); % Plots only the shapes with no stroke
                    if aSlice ==n
                        % store patch for legend
                        P.allObj = [P.allObj pp.(ss)];
                    end
                    uistack(pp.(ss),'top')
                    %     if strcmp(ss,'s1') % Adding a line only for the model
                    %         tt = plot(X.(ss)(2:end-1,aSlice),Y.(ss)(2:end-1,aSlice),'-','color',col*3/4,'linewidth',0.75); % Plots only the shapes with no stroke
                    %         uistack(tt,'top')
                    %     end
                    FaceAlpha = storeFaceAlpha;
                case 'line' 
                    col = cmap(thisRowCol,:);
                    pp.(ss) = plot(X.(ss)(:,aSlice),Y.(ss)(:,aSlice),'-','color',col,'linewidth',2); % Plots only the shapes with no stroke
            end
        end
    end
    if nS ==1
        % Y-Axis labels
        ylabel(['$\mathrm{\mathbf{Time \; (d)}}$'],'interpreter','latex')
        yticks(locationHorizLines(1,:))
        xlabel('\textbf{Microcolony size ($\mathbf{\mu} \mathrm{\mathbf{m}}\mathbf{^3}$)}','Interpreter','latex') 
    end 
    if VLines ~= 0
        locationHorizLines_inv = fliplr(locationHorizLines);
        plotStatistics(x,mini,maxi,data_s,showSlices_s,aSlice,ss,locationHorizLines_inv,VLinesWidth,[])
    end
 end
end

%% Finalizing appearance
h = gca();
interPlotHeight = diff(h.YTick); 
switch IC_at
 case 'SS_computed'
    heightIQbars = interPlotHeight(end);
 case 'startupScenario'
    heightIQbars = maxHeightExpData;
end
                    
%% adding vertical lines at the grid points
for zz = 1:length(x)
 rr = plot([x(zz) x(zz)],h.YLim,'-','color',[0.8 0.8 0.8],'linewidth',0.3);
 uistack(rr,'bottom');
end 
%% adding horizLines lines 
for zz = 1:size(locationHorizLines,2)
 rr = plot(h.XLim,[locationHorizLines(1,zz) locationHorizLines(1,zz)],'-','color',[0.8 0.8 0.8],'linewidth',0.3);
 uistack(rr,'bottom');
end 

%%
objj = fig1.findobj;
objj(2).CLim = [0.3 0.96];
% set(gca, 'ylim',[0 5.5])
set(gca,'FontSize',17)
set(gca,'LineWidth',1)
set(gca,'TickDir','out'); % The only other option is 'in'
set(ax,'TickLength',[0.005 0.005])
xtickpos = get(gca, 'xtick');
xtLabels = get(gca, 'XTickLabel');
% xtickpos = xtickpos+mini;
xticks(xtickpos)
xticklabels(xtLabels);

ytickpos = get(gca, 'ytick');
ylbls = string(yTickLabels);
if ~reverse 
 ylbls = flip(ylbls);
end
switch IC_at 
 case 'startupScenario'
    % removing some labels otherwise too dense
    if ~reverse 
        ylbls([end-1,end-2]) = [];
    else
        ylbls([2,3]) = [];
    end
    ytickpos([end-1,end-2]) = [];
end
yticks(ytickpos)
yticklabels(ylbls)
yyaxis right
yticks([])
yticklabels([])
ylabel('NDF')
yyaxis left
h.XLim = [0 maxi];
switch overlapMethod
 case 'variable'
    h.YLim = [0 h.YTick(end)+ heightIQbars*3];
 case 'constant'
    h.YLim = [0 h.YTick(end)+ heightIQbars];
end

hold off

fff = gca;
fff.YAxis(2).Color = hex2rgb('7E2F8E');

if iT ==noTrials 
% annotation('textbox', [0.51, 0.89, 0.24, 0.03], 'String','\textbf{SRT }$\mathbf{ = 20}$\textbf{ (bottom)}$\mathbf{,15,10,5}$\textbf{ (top) \textit{d}}','Interpreter','latex','FitBoxToText','on','BackgroundColor','w','EdgeColor',[0.8 0.8 0.8],'LineWidth',0.5,'Margin',2,'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','FontSize',14);
 %annotation('textbox', [0.58, 0.88, 0.3, 0.03], 'String','SRT = 20 (bottom),15,10,5 (top)','FitBoxToText','on','BackgroundColor','w','EdgeColor','k','LineWidth',0.5,'Margin',2,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
 %delete(findall(gcf,'type','annotation')) 
 P.allObjNames(end) = strrep(P.allObjNames(end),'d;','d.'); 
 
 [hleg,hobj] = legend(P.allObj,P.allObjNames,'Interpreter','latex','location','northeast','Orientation','horizontal');   
 for p = (noTrials+1)+1:2*(noTrials+1)
    hobj(p).Vertices(1:2,1) = hobj(p).Vertices(1:2,1)+0.05;
    hobj(p).Vertices(:,2) = hobj(p).Vertices(:,2)+0.02;
 end
 aa_ll = 0;
 
 for p = (noTrials+1):-1:1
    aa_ll = aa_ll+0.04;
    hobj(p).Position(1) = hobj(p).Position(1)+aa_ll;
 end
 aa_ll = 0;
 for p = 2*(noTrials+1):-1:(noTrials+1)+1
    aa_ll = aa_ll+0.04;
    hobj(p).Vertices(:,1) = hobj(p).Vertices(:,1)+aa_ll;
 end
 title(hleg,'$\mathrm{\mathbf{SRT = }}$','Interpreter','latex','FontSize',18);
 hleg.Title.Visible = 'on';
 hleg.Title.NodeChildren.Position = [0.06 0.52 0.6];
 legend('boxoff')
 
 set(gcf, 'Renderer', 'painters')

 saveas(fig1,['./outputs/figures/joyPlot/joy_lambda_' labLambda improvStr cycleModStr '_lowSRT.fig'])
 saveas(fig1,['./outputs/figures/joyPlot/joy_lambda_' labLambda improvStr cycleModStr '_lowSRT.pdf']) 
end

joyFigs.fig1 = fig1;
return

function [] = plotStatistics(x,mini,maxi,data_s,showSlices_s,aSlice,ss,locationHorizLines_inv,VLinesWidth,Ylim)
switch ss
 case 's1'
    % model
    colModel = 'k';
    if ismember(aSlice,showSlices_s.(ss))
    % mean 
        yv = interp1(x,data_s.s1.model(:,aSlice),data_s.s1.mean_model(aSlice)); % from top to bottom
        XV = [data_s.s1.mean_model(aSlice) data_s.s1.mean_model(aSlice)];
        YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
        if isempty(XV(XV<mini|XV>maxi))
            plot(XV,YV,'-','color',colModel,'LineWidth',VLinesWidth);
        end
    % q1 
        yv = interp1(x,data_s.s1.model(:,aSlice),data_s.s1.q1_model(aSlice)); % from top to bottom
        XV = [data_s.s1.q1_model(aSlice) data_s.s1.q1_model(aSlice)];
        % yv = heightIQbars;
        YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
        if isempty(XV(XV<mini|XV>maxi))
            plot(XV,YV,':','color',colModel,'LineWidth',VLinesWidth*5/3);
        end

    % q3 
        yv = interp1(x,data_s.s1.model(:,aSlice),data_s.s1.q3_model(aSlice)); % from top to bottom
        XV = [data_s.s1.q3_model(aSlice) data_s.s1.q3_model(aSlice)];
        % yv = heightIQbars;
        YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
        if isempty(XV(XV<mini|XV>maxi))
            plot(XV,YV,':','color',colModel,'LineWidth',VLinesWidth*5/3);
        end
    end
 case 's2'
    colData = 'r';
    if ismember(aSlice,showSlices_s.(ss))
        % data
        % mean 
            yv = interp1(x,data_s.s2.data(:,aSlice),data_s.s2.mean_data(aSlice)); % from top to bottom
            XV = [data_s.s2.mean_data(aSlice) data_s.s2.mean_data(aSlice)];
            YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
            if ~isempty(Ylim)
                YV = Ylim;
            end
            if isempty(XV(XV<mini|XV>maxi))
                plot(XV,YV,'-','color',colData,'LineWidth',VLinesWidth);
            end
 % 
        % q1 
            yv = interp1(x,data_s.s2.data(:,aSlice),data_s.s2.q1_data(aSlice)); % from top to bottom
            XV = [data_s.s2.q1_data(aSlice) data_s.s2.q1_data(aSlice)];
            % yv = heightIQbars;
            YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
            if ~isempty(Ylim)
                YV = Ylim;
            end
            if isempty(XV(XV<mini|XV>maxi))
                plot(XV,YV,':','color',colData,'LineWidth',VLinesWidth*5/3);
            end

        % q3 
            yv = interp1(x,data_s.s2.data(:,aSlice),data_s.s2.q3_data(aSlice)); % from top to bottom
            XV = [data_s.s2.q3_data(aSlice) data_s.s2.q3_data(aSlice)];
            % yv = heightIQbars;
            YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
            if ~isempty(Ylim)
                YV = Ylim;
            end
            if isempty(XV(XV<mini|XV>maxi))
                plot(XV,YV,':','color',colData,'LineWidth',VLinesWidth*5/3)
            end      
    end
    
  case 's5'
    if ismember(aSlice,showSlices_s.(ss))
     % mean 
        yv = interp1(x,data_s.(ss).model(:,aSlice),data_s.(ss).mean_model(aSlice)); % from top to bottom
        XV = [data_s.(ss).mean_model(aSlice) data_s.(ss).mean_model(aSlice)];
        YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
        if isempty(XV(XV<mini|XV>maxi))
            plot(XV,YV,'--','color',[91, 207, 244] / 255,'LineWidth',VLinesWidth);
        end
    % q1 
        yv = interp1(x,data_s.(ss).model(:,aSlice),data_s.(ss).q1_model(aSlice)); % from top to bottom
        XV = [data_s.(ss).q1_model(aSlice) data_s.(ss).q1_model(aSlice)];
        % yv = heightIQbars;
        YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
        if isempty(XV(XV<mini|XV>maxi))
            plot(XV,YV,':','color',[91, 207, 244] / 255,'LineWidth',VLinesWidth*5/3);
        end

    % q3 
        yv = interp1(x,data_s.(ss).model(:,aSlice),data_s.(ss).q3_model(aSlice)); % from top to bottom
        XV = [data_s.(ss).q3_model(aSlice) data_s.(ss).q3_model(aSlice)];
        % yv = heightIQbars;
        YV = [locationHorizLines_inv(1,aSlice);locationHorizLines_inv(1,aSlice) + yv];
        if isempty(XV(XV<mini|XV>maxi))
            plot(XV,YV,':','color',[91, 207, 244] / 255,'LineWidth',VLinesWidth*5/3);
        end
    end
end
return


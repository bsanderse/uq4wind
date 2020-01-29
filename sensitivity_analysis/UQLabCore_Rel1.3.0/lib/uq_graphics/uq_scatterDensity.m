function [h,Axs,bigAx] = uq_scatterDensity(varargin)
%UQ_SCATTERDENSITY creates a scatter plot matrix of multivariate data.
%    
%   UQ_SCATTERDENSITY(X) creates a scatter plot matrix from each pairs of
%   the columns in X. If X is N-by-M matrix, UQ_SCATTERDENSITY produces a
%   lower triangular matrix of M-by-M plots, in which histograms of each
%   column of X are plotted in the main diagonal and scatter plots of each
%   pair of columns of X are plotted as the off-diagonal elements. 
%   
%   UQ_SCATTERDENSITY(..., hist_NAME, VALUE, scatter_NAME, VALUE) modifies
%   the scatter plot matrix using additional NAME/VALUE pair arguments. The
%   following convention applies when providing pairs to the corresponding
%   elements of the matrix of plots:
%       - If NAME has the prefix 'hist_', then VALUE is passed to the
%         function <a href="matlab:help uq_histogram"
%         >uq_histogram</a> to create the histograms in the main diagonal.
%       - If NAME has the prefix 'scatter_', then VALUE is passed instead
%         to the <a href="matlab:help uq_plot"
%         >uq_plot</a> function that creates the off-diagonal scatter plots.
%
%   UQ_SCATTERDENSITY(..., NAME, VALUE) modifies the overall plot using 
%   NAME/VALUE pair arguments of which the NAMEs have neither the prefix
%   'hist_' nor 'scatter_'. Available options are:
%   
%      Name               VALUE
%      'points'           A set of N points that are plotted together with
%                         the X points. 
%                         - Numerical N-by-M
%
%      'labels'           A cell array of labels for subplots
%                         - Cell Array 1-by-M
%
%      'color'            A vector specifying the color of the plotted
%                         elements
%                         - Numerical 1-by-3
% 
%      'title'            Title string to add to plot
%
%   UQ_SCATTERDENSITY(AX,...) creates the plots into the Axes object AX
%   instead of the current axes.
%
%   H = UQ_SCATTERDENSITY(...) returns the graphics objects of the plot.
%   The main diagonal elements (histograms) are represented by Bar objects
%   and the off-diagonal elements are represented by Line objects. Use the 
%   elements in H to access and modify the properties of specific plots. In
%   MATLAB R2014a or older, the function returns handles to barseries
%   objects (instead of Bar) and lineseries objects (instead of Line).
%
%   [H,AXS] = UQ_SCATTERDENSITY(...) additionally returns the Axes objects
%   (or handles) of all the plots inside the scatter plot matrix.
%
%   [H,AXS,BIGAX] = UQ_SCATTERDENSITY(...) additionally returns the parent
%   axes of the scatter plot matrix.
%
%   See also UQ_HISTOGRAM, UQ_PLOT.

%% Verify inputs
if nargin == 0
    error('Not enough input arguments.')
end

%% Default plot properties
DEFAULTS.color = uq_colorOrder(1);

% Histogram defaults
hist_DEFAULTS.EdgeColor = 'none';  % No edge color
hist_DEFAULTS.FaceAlpha = 1;       % No transparancy

% scatter plot defaults
scatter_DEFAULTS.MarkerSize = 5;
scatter_DEFAULTS.Marker = '.';
scatter_DEFAULTS.LineStyle = 'none';

% OS-specific properties
if ispc % windows
    
elseif isunix||ismac % linux
    
end

%% Get the axes to plot in
ParentAx = uq_getPlotAxes(varargin{:});
% Set the axes to be invisible
axis(ParentAx,'off')

%% Parse inputs, split NAME/VALUE pairs
isAxes = uq_isAxes(varargin{1});

if isAxes
    if isnumeric(varargin{2})
        X = varargin{2};
        [histArg,scatterArg,uqArg] = uq_splitScatterArgs(varargin(3:end));
    else
        error('Data must be provided.')
    end
elseif isnumeric(varargin{1})
    X = varargin{1};
    [histArg,scatterArg,uqArg] = uq_splitScatterArgs(varargin(2:end));
else
    error('Either data or axes object must be provided as the first argument.')
end

%% Parse NAME/VALUE pairs other than histogram and scatter
parse_keys = {'points','labels','color','title'};
parse_types = {'p','p','p','p'};
% make NAME lower case
uqArg(1:2:end) = lower(uqArg(1:2:end));
[uq_cline,~] = uq_simple_parser(uqArg, parse_keys, parse_types);

% 'points' option
if ~strcmp(uq_cline{1}, 'false')
    plotPoints_flag = true;
    plotPoints = uq_cline{1};
    % Check dimensions
    if ~(size(X,2) == size(plotPoints,2))
        error('Additional plot points do not have a compatible size')
    end
else
    plotPoints_flag = false;
end

% 'labels' option
if ~strcmp(uq_cline{2}, 'false')
    plotLabels_flag = true;
    plotLabels = uq_cline{2};
    % Check dimensions
    if ~(size(X,2) == length(plotLabels))
        error('Label length does not match the provided points dimension')
    end
else
    plotLabels_flag = false;
end

% 'color' option
if ~strcmp(uq_cline{3}, 'false')
    plotColor = uq_cline{3};
else
    plotColor = DEFAULTS.color;
end
% Use specified 'color' in histogram and scatter plot if they are not
% specified.
if isempty(histArg) || ~any(any(strcmpi(histArg,'facecolor')))
    histArg(:,end+1) = {'FaceColor'; plotColor};
end
if isempty(scatterArg) || ~any(any(strcmpi(scatterArg,'color')))
    scatterArg(:,end+1) = {'MarkerEdgeColor'; plotColor};
end

% 'title' option
if ~strcmp(uq_cline{4}, 'false')
    titleString = uq_cline{4};
    TitleProp = get(ParentAx,'Title');
    set(TitleProp,...
        'String', titleString,...
        'FontSize', 16,...
        'Interpreter', 'LaTeX');
end

%% Get the number of plots
nDim = size(X,2); % Number of plots

%% Compute the font size scaling factor
% Scale font size based on number of plots linearly between 1 (up to 2 
% subplots) and 0.2 (more than 10 subplots)
lowerPlots = 2;
upperPlots = 10;
lowerScale = 1;
upperScale = 0.3;
fontScale = interp1(...
    [1,lowerPlots,upperPlots], [lowerScale,lowerScale,upperScale],...
    nDim, 'linear', upperScale);  

%% Start plotting

% Get smallest and largest sample values in each dimension
% to set axis properly
if plotPoints_flag
    minSample = min([X;plotPoints]);
    maxSample = max([X;plotPoints]);
else
    minSample = min(X);
    maxSample = max(X);
end
rangeSample = maxSample - minSample;

% Compute plot limits
plotPadding = 0.05;
plotLimits = [minSample - plotPadding*rangeSample;...
              maxSample + plotPadding*rangeSample;]; 

% Loop over dimensions to plot
for ii = 1:nDim
    for jj = 1:nDim
        % Switch between 'uq_histogram' and 'uq_plot'
        if ii == jj %histogram
            % Create axes
            currSub = uq_subplot(nDim,nDim,(ii-1)*nDim+jj);
            % Plot histogram 
            currObj = uq_histogram(currSub, X(:,ii), histArg{:});
            % Add additional formatting
            uq_formatGraphObj(currObj, histArg, hist_DEFAULTS)
            currSubFontSize = get(currSub,'FontSize');
            axes_DEFAULTS.FontSize = fontScale*currSubFontSize;
            uq_formatGraphObj(currSub, [], axes_DEFAULTS)

            % Display point estimate if requested
            if plotPoints_flag
                hold on
                uq_plot(...
                    currSub,...
                    [plotPoints(:,ii) plotPoints(:,ii)],...
                    get(gca,'ylim'), 'r')
                hold off
            end
        elseif jj < ii
            % Create axes
            currSub = uq_subplot(nDim,nDim,(ii-1)*nDim+jj);
            % Plot scatter
            currObj = uq_plot(currSub, X(:,jj), X(:,ii), scatterArg{:});
            % Add default formatting
            uq_formatGraphObj(currObj, scatterArg, scatter_DEFAULTS)
            currSubFontSize = get(currSub,'FontSize');
            axes_DEFAULTS.FontSize = fontScale*currSubFontSize;
            uq_formatGraphObj(currSub, [], axes_DEFAULTS)
            % Display point estimate also
            if plotPoints_flag
                hold on
                uq_plot(...
                    plotPoints(:,jj), plotPoints(:,ii), 'r+',...
                    'MarkerSize', 8, 'LineStyle', 'none')
                hold off
            end           
        end
        
        % Take care of axis labels and ticks
        NumTicks = 2;
        xTicks = linspace(minSample(jj),maxSample(jj),NumTicks);
        yTicks = linspace(minSample(ii),maxSample(ii),NumTicks);
        if ii == jj
            % Diagonal elements
            set(currSub,...
                'XTick', xTicks,...
                'YTick', [],...
                'XTickLabel', [],...
                'YTickLabel', [],...
                'XLim', plotLimits(:,ii)');
            if jj == 1
                % First element of the diagonal
                if plotLabels_flag
                    ylabel(currSub,plotLabels{ii})
                end
            end
            if jj == nDim
                % Last element of the diagonal
                if plotLabels_flag
                    xlabel(currSub,plotLabels{ii})
                end
                set(currSub, 'XTickLabel', xTicks)
                % NOTE: for backward compatibility with R2014a
                if which('xtickformat')
                    xtickformat('%.2f')
                end
                % NOTE: for backward compatibility with R2014a
                if which('xtickangle')
                    xtickangle(45)
                end
            end
        elseif jj < ii
            % Off-Diagonal elements
            set(currSub,...
                'YTick',yTicks,...
                'XTick',xTicks,...
                'XTickLabel', [],...
                'YTickLabel', [],...
                'XLim', plotLimits(:,jj)',...
                'YLim', plotLimits(:,ii)');
            if jj == 1
                % First column (axes and tick labels on the left)
                if plotLabels_flag
                    ylabel(currSub,plotLabels{ii})
                end
                set(currSub,'YTickLabel',yTicks)
                % NOTE: for backward compatibility with R2014a
                if which('ytickformat')
                    ytickformat('%.2f')
                end
                % NOTE: for backward compatibility with R2014a
                if which('ytickangle')
                    ytickangle(45)
                end
            end
            if ii == nDim
                % Last row (axes and tick labels on the bottom)
                if plotLabels_flag
                    xlabel(currSub,plotLabels{jj})
                end
                set(currSub,'XTickLabel',xTicks)
                % NOTE: for backward compatibility with R2014a
                if which('xtickformat')
                    xtickformat('%.2f')
                end
                % NOTE: for backward compatibility with R2014a
                if which('xtickangle')
                    xtickangle(45)
                end
            end
        end
      
        % Store axes and objects
        AX(ii,jj) = currSub;
        S(ii,jj) = currObj;
    end
end

if plotPoints_flag
    % Add callback for figure change to scale point estimate line
    try
        % add callback to figure
        set(currFig,'SizeChangedFcn',{@resizeui, pointLine, AX});
        notify(currFig,'SizeChanged')
    catch
        % not supported
    end
end

%% Return the output if requested
if nargout > 0
    h = S;
    Axs = AX;
    bigAx = ParentAx;
end

end

%% ------------------------------------------------------------------------
function resizeui(hObject ,event, LineContainer, AxesContainer)
% loop over lines and set y limit dynamically if the figure is resized

for ii = 1:length(LineContainer)
    yLimCurr = AxesContainer(ii,ii).YLim;
    set(LineContainer(ii), 'YData', yLimCurr);
end

end

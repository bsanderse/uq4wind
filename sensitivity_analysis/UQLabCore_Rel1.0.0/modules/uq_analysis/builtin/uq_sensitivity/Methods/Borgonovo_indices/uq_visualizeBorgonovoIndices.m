function uq_visualizeBorgonovoIndices(Results, Outputs, varargin)
% UQ_VISUALIZEBORGONOVOINDICES(RESULTS,OUTIDX,VARARGIN): graphical
%     representation of the Beorgonovo indices in the results structure
%     RESULTS.
%
% See also: UQ_DISPLAY_UQ_SENSITIVITY

%% Universal options for the Sobol' Indices Plots
GlobalOptions.NIndices = 10; % max. number of Indices to plot 
GlobalOptions.NIndicesHigh = 5; % max number of indices to plot (higher orders)
GlobalOptions.YTicks = 0:0.2:1; % Ticks of y-axis in plots
GlobalOptions.AxisFontSize = 14; 
GlobalOptions.NumOfColors = 6;
GlobalOptions.IdxSorting = 'descend'; % How to sort high-order indices
GlobalOptions.BoxBorder = 1.5; % Thickness of box around plot

%% parse varargin options
% initialization
pie_flag = false;
hist_flag = true;
joint_flag = false;

if nargin > 2
    parse_keys = {'pie', 'hist','Joint PDF'};
    parse_types = {'f', 'f','f'};
    [uq_cline, varargin] = uq_simple_parser(varargin, parse_keys, parse_types);
    
    % 'coefficients' option additionally prints the coefficients
    if strcmp(uq_cline{1}, 'true')
        pie_flag = true;
        hist_flag = false;
    end
    
    % 'tolerance' option sets the default tolerance to plot coefficients
    if strcmp(uq_cline{2}, 'true')
        hist_flag = true;
    end
    
    % This will result in a visualization of the estimate of the joint
    % distribution used:
    if strcmp(uq_cline{3}, 'true')
        joint_flag = true;
        joint_to_vis = varargin{1};
    end

end



%% Parse Input Arguments
if ~exist('Outputs', 'var') || isempty(Outputs)
    Outputs = 1;
end

Tag = '';
%Bootstrap bounds currently not supported:
Bootstrap = 0;

if nargin < 5 % Maybe some consistency checks...
    MyColors = uq_cmap(GlobalOptions.NumOfColors);
else
    MyColors = varargin{4} ;
end
%% Produce plot(s)
% Check how many outputs are there:
NOuts = length(Outputs);
OutputTag = '';
for oo = Outputs
    % Get the number of indices
    NumOfIndices = length(Results.Delta(:,oo));
    if NOuts > 1
        OutputTag = sprintf(', output #%d', oo);
    end
    
    % --- 1st order INDICES --- %
    if hist_flag && ~joint_flag
        % create the figure ("1st Order" deleted since only first order
        % available)
        uq_figure('Name',...
            sprintf('%sBorgonovo indices %s', OutputTag, Tag),...
            'Position', [50 50 500 400]);
        hold on
        % create the bar plot
        x = 1:NumOfIndices;
        uq_bar(x, Results.Delta(:,oo),...
            'FaceColor', MyColors(1,:), 'EdgeColor', 'none');
        % Fix labels, title and axis properties of the plot
        uq_setInterpreters(gca)
        customize_Borgonovo_plot(gca, oo, Results, 'delta', GlobalOptions)
        set(gca, 'LineWidth', GlobalOptions.BoxBorder)
    end
    
    if pie_flag
        % create the figure
        uq_figure('name', sprintf('%sBorgonovo indices (pie)%s', Tag,OutputTag, 'Position', [50 50 500 400]));
        % create the bar plot
        x = zeros(size(Results.Delta(:,oo)));
        pp = pie3(Results.Delta(:,oo), x, cellstr(Results.VariableNames'));
        set(findobj(pp, 'type', 'surface'), 'edgecolor', 'k')
        set(findobj(pp, 'type', 'patch'), 'edgecolor', 'k')
        material metal
        camlight
        rotate3d on
        set(gca, 'LineWidth', GlobalOptions.BoxBorder)
    end
    
    if joint_flag
        % Creates a visualization of the joint pdf used as it was used in 
        % the computation:
        show_conditionalPDFestimates(Results,oo,joint_to_vis,GlobalOptions)
    end

end

end

function customize_Borgonovo_plot(ax, oo, Results, IndType,...
    GlobalOptions, current_order, idx)
    
switch lower(IndType)
    case 'delta'
       NumOfIndices = length(Results.Delta(:,oo));
       title('Borgonovo indices')
       ylabel('$\mathrm{\delta_i}$') % "Order 1" removed
       customize_axes(ax, oo, Results, IndType, GlobalOptions)
       ylim([0 1])
end


xlim([0, NumOfIndices + 1]);
box on;
end

function customize_axes(ax, oo, Results, IndType, GlobalOptions,...
    current_order, idx)
% Helper function for setting the xticks of an axis object used
% for plotting of Sobol' indices
NIndices = GlobalOptions.NIndices;
YTicks = GlobalOptions.YTicks;
fs = GlobalOptions.AxisFontSize;

switch lower(IndType)
    case 'delta'
        NumOfIndices = length(Results.Delta(:,oo));
        XTicks = ceil(linspace(1,NumOfIndices, min(NumOfIndices,NIndices)));
        TickNames = Results.VariableNames(XTicks);
end

set(ax, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
    'FontSize', fs)

end

function show_conditionalPDFestimates(Results, oo, idx, GlobalOptions)

if isnumeric(idx) && ...
        max(idx)>length(Results(oo).Delta)
    error('The requested variable index does not exist in the model!');
end

if strcmpi(idx,'all')
    idx = 1:length(Results.Delta(:,oo));
end

for ii = idx
    uq_figure('Name',...
        sprintf('Joint Pdf Y|%s',Results.VariableNames{ii}),...
        'Position', [50 50 500 400])
    pdfplot = pcolor(Results.JointPDF{ii,oo}');
    hold on
    shading flat
    set(pdfplot, 'EdgeColor', 'none')
    box on
%     set(pdfplot,'EdgeAlpha',0)
    % Set labels
    s = sprintf('%s',Results.VariableNames{ii});
    set(gca, 'FontSize', 14)
    uq_setInterpreters(gca)
    xlabel(s)
    ylabel('Y')
    set(gca, 'LineWidth', 2)
    hold off
%     title(sprintf('f_{Y,X_%i}(x_%i,y)',ii,ii));
end
%figure;
%surf(ccdf);
end
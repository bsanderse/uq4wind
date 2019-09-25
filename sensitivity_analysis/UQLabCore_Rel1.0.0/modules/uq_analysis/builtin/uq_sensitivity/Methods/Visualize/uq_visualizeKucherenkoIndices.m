function uq_visualizeKucherenkoIndices(Results, Outputs, varargin)
% UQ_KUCHERENKOINDICES(RESULTS,OUTIDX,VARARGIN): graphical
%     representation of the Kucherenko indices in the results structure
%     RESULTS.
%
% See also: UQ_DISPLAY_UQ_SENSITIVITY

%% Universal options for the ANCOVA Indices Plots
GlobalOptions.NIndices = 10; % max. number of Indices to plot (total and 1st order)
% GlobalOptions.NIndicesHigh = 5; % max number of indices to plot (higher orders)
GlobalOptions.YTicks = 0:0.2:10; % Ticks of y-axis in plots
GlobalOptions.AxisFontSize = 14; 
GlobalOptions.NumOfColors = 6;
% GlobalOptions.IdxSorting = 'descend'; % How to sort high-order indices

%% Parse Input Arguments
if ~exist('Outputs', 'var') || isempty(Outputs)
    Outputs = 1;
end

Tag = '';

if nargin < 4 % Maybe some consistency checks...
    MyColors = uq_cmap(GlobalOptions.NumOfColors);
else
    MyColors = varargin{2} ;
end
%% Produce plot(s)
% Check how many outputs are there:
NOuts = length(Outputs);
OutputTag = '';

% Grab earlier specifications
NIndices = GlobalOptions.NIndices;
YTicks = GlobalOptions.YTicks;
fs = GlobalOptions.AxisFontSize ;

for oo = Outputs
    % Get the number of indices
    NumOfIndices = length(Results.Total(:,oo));
    XTicks = ceil(linspace(1,NumOfIndices, min(NumOfIndices,NIndices)));
    TickNames = Results.VariableNames(XTicks);
    if NOuts > 1
        OutputTag = sprintf(', output #%d', oo);
    end
    
    % --- TOTAL INDICES --- %
    % Create the figure
    uq_figure('name',...
        sprintf('%sTotal Kucherenko indices%s', Tag, OutputTag),...
        'Position', [50 50 500 400])
    hold on
    % Create the bar plot
    x = 1:NumOfIndices;
    uq_bar(x, Results.Total(:,oo),...
        'FaceColor', MyColors(1, :), 'EdgeColor', 'none')
    
    % Fix title, labels, and axis properties of the plot
    uq_setInterpreters(gca)
    title('Total Kucherenko indices')
    ylabel('$\mathrm{S_i^{T}}$')
    set(gca, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
        'FontSize', fs, 'YTickMode', 'auto','Linewidth',2)
    
    % --- FIRST ORDER INDICES --- %
    % Create the figure
    uq_figure('name',...
        sprintf('%sFirst-order Kucherenko indices%s', Tag, OutputTag),...
        'Position', [50 50 500 400])
    hold on
    % Create the bar plot
    x = 1:NumOfIndices;
    uq_bar(x, Results.FirstOrder(:,oo),...
        'FaceColor', MyColors(2, :), 'EdgeColor', 'none')

    % Fix title, labels, axis properties of the plot
    uq_setInterpreters(gca)
    title('First-order Kucherenko indices')
    ylabel('$\mathrm{S_i}$');
    set(gca, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
        'FontSize', fs, 'YTickMode', 'auto', 'Linewidth', 2)
    
end

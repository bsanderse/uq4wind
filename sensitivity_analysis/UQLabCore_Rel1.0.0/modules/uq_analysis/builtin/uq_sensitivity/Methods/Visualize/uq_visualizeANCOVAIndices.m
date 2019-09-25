function uq_visualizeANCOVAIndices(Results, Outputs, varargin)
% UQ_VISUALIZEANCOVAINDICES(RESULTS,OUTIDX,VARARGIN): graphical
%     representation of the ANCOVA indices in the results structure
%     RESULTS.
%
% See also: UQ_DISPLAY_UQ_SENSITIVITY

%% Universal options for the ANCOVA Indices Plots
GlobalOptions.NIndices = 10; % max. number of Indices to plot (total and 1st order)
GlobalOptions.YTicks = 0:0.2:10; % Ticks of y-axis in plots
GlobalOptions.AxisFontSize = 14; 
GlobalOptions.NumOfColors = 6;

% %% parse varargin options
% % initialization
% pie_flag = false;
% hist_flag = true;
% 
% if nargin > 2
%     parse_keys = {'pie', 'hist'};
%     parse_types = {'f', 'f'};
%     [uq_cline, varargin] = uq_simple_parser(varargin, parse_keys, parse_types);
%     % 'coefficients' option additionally prints the coefficients
%     if strcmp(uq_cline{1}, 'true')
%         pie_flag = true;
%         hist_flag = false;
%     end
%     
%     % 'tolerance' option sets the default tolerance to plot coefficients
%     if strcmp(uq_cline{2}, 'true')
%         hist_flag = true;
%     end
% 
% end

%% Parse Input Arguments
if ~exist('Outputs', 'var') || isempty(Outputs)
    Outputs = 1;
end
% 
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
    NumOfIndices = length(Results.FirstOrder(:,oo));
    XTicks = ceil(linspace(1,NumOfIndices, min(NumOfIndices,NIndices)));
    TickNames = Results.VariableNames(XTicks);
    
    if NOuts > 1
        OutputTag = sprintf(', output #%d', oo);
    end
    
    % Set the y-axis limits and make sure to have two different values
    round_limit = 0.15;
    ylimitunc = ([floor(min(Results.Uncorrelated(:,oo))/round_limit)*round_limit,...
        ceil(max(Results.Uncorrelated(:,oo))/round_limit)*round_limit]);
    if ylimitunc(1)==ylimitunc(2)
        ylimitunc(2) = ylimitunc(1)+0.3;
    end
    ylimitint = ([floor(min(Results.Interactive(:,oo))/round_limit)*round_limit,...
        ceil(max(Results.Interactive(:,oo))/round_limit)*round_limit]);
    if ylimitint(1)==ylimitint(2)
        ylimitint(2) = ylimitint(1)+0.3;
    end
    ylimitcor = ([floor(min(Results.Correlated(:,oo))/round_limit)*round_limit,...
        ceil(max(Results.Correlated(:,oo))/round_limit)*round_limit]);
    if ylimitcor(1)==ylimitcor(2)
        ylimitcor(2) = ylimitcor(1)+0.3;
    end
    ylimitfis = ([floor(min(Results.FirstOrder(:,oo))/round_limit)*round_limit,...
        ceil(max(Results.FirstOrder(:,oo))/round_limit)*round_limit]);
    if ylimitfis(1)==ylimitfis(2)
        ylimitfis(2) = ylimitfis(1)+0.3;
    end
    
    % --- UNCORRELATED INDICES --- %
    % Create the figure
    uq_figure('name',...
        sprintf('%sUncorrelated ANCOVA indices%s', Tag, OutputTag),...
        'Position', [50 50 500 400])
    hold on
    % Create the bar plot
    x = 1:NumOfIndices;
    uq_bar(x, Results.Uncorrelated(:,oo),...
        'FaceColor', MyColors(1,:), 'EdgeColor', 'none')
    ylim(ylimitunc)

    % Fix labels, title and axis properties of the plot
    uq_setInterpreters(gca)
    title('Uncorrelated ANCOVA indices')
    ylabel('$\mathrm{S_i^{U}}$')
    set(gca, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
        'FontSize', fs, 'YTickMode', 'auto', 'LineWidth', 2)

    % --- INTERACTIVE indices --- %
    % Create the figure
    uq_figure('name',...
        sprintf('%sInteractive ANCOVA indices%s', Tag, OutputTag),...
        'Position', [50 50 500 400])
    hold on
    % Create the bar plot
    x = 1:NumOfIndices;
    uq_bar(x, Results.Interactive(:,oo),...
        'FaceColor', MyColors(2,:), 'EdgeColor', 'none');
    ylim(ylimitint) 
    
    % Fix labels, title and axis properties of the plot
    uq_setInterpreters(gca)
    title('Interactive ANCOVA indices')
    ylabel('$\mathrm{S_i^{I}}$')
    set(gca, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
        'FontSize', fs, 'YTickMode', 'auto', 'Linewidth', 2);
       
    % --- CORRELATED INDICES --- %
    % Create the figure
    uq_figure('name',...
        sprintf('%sCorrelated ANCOVA indices%s', Tag, OutputTag),...
        'Position', [50 50 500 400])
    hold on
    % Create the bar plot
    x = 1:NumOfIndices;
    uq_bar(x, Results.Correlated(:,oo),...
        'FaceColor', MyColors(3,:), 'EdgeColor', 'none')
    ylim(ylimitcor)

    % Fix titles, labels, and axis properties of the plot
    uq_setInterpreters(gca)
    title('Correlated ANCOVA indices')
    ylabel('$\mathrm{S_i^{C}}$')
    set(gca, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
        'FontSize', fs, 'YTickMode', 'auto', 'Linewidth', 2)
    
    % --- SUMMED UP FIRST ORDER INDICES --- %
    % Create the figure
    uq_figure('name',...
        sprintf('%sFirst order ANCOVA indices%s', Tag, OutputTag),...
        'Position', [50 50 500 400])
    hold on
    % Create the bar plot
    x = 1:NumOfIndices;
    uq_bar(x, Results.FirstOrder(:,oo),...
        'FaceColor', MyColors(4,:), 'EdgeColor', 'none')
    ylim(ylimitfis)
    
    % Fix title, labels, and axis properties of the plot
    uq_setInterpreters(gca)
    title('First-order ANCOVA indices')
    ylabel('$\mathrm{S_i}$')
    set(gca, 'XTick', XTicks, 'XTickLabel', TickNames, 'YTick', YTicks,...
        'FontSize', fs, 'YTickMode', 'auto','Linewidth',2)
    
end
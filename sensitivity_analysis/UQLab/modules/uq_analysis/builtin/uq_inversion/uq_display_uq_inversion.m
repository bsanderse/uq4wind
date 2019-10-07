function uq_display_uq_inversion(module, varargin)
% UQ_DISPLAY_UQ_INVERSION graphically displays the results of an inverse
%    analysis carried out with the Bayesian inversion module of UQLab.
%
%    UQ_DISPLAY_UQ_INVERSION(MODULE, NAME, VALUE) allows to choose
%    more advanced plot functions by specifying Name/Value pairs:
%
%       Name               VALUE
%       'acceptance'       Plots the acceptance rate per chain 
%                          - Logical
%                          default : false
%       'scatterplot'      Plots a multi dimensional parameter scatter plot
%                          of the parameters in the Results field of MODULE 
%                          - Integer or 'all'
%                          default : 'all'
%       'predDist'         Plots the predictive distributions if available
%                          from post processing (uq_postProcessInversion)
%                          - Logical
%                          default : true (if available)
%       'meanConvergence'  Plots the convergence of the marginal mean for
%                          the specified parameter averaged over all chains 
%                          - Integer or 'all'
%                          default : 'none'
%       'trace'            Plots the trace plot of the marginal sample
%                          points for the specified parameter 
%                          - Integer or 'all'
%                          default : 'none'
%                          
% See also: UQ_PRINT_UQ_INVERSION, UQ_POSTPROCESSINVERSION

%% CONSISTENCY CHECKS
if ~strcmp(module.Type, 'uq_inversion')
   fprintf('uq_display_uq_inversion only operates on objects of type ''Inversion''') 
end

%check if MCMC Solver
if or(~strcmp(module.Internal.Solver.Type, 'MCMC'),isempty(module.Results))
    error('No results to display')
end

% switch if custom likelihood
if module.Internal.customLikeli
    CUSTOM_LIKELI = true;
else
    CUSTOM_LIKELI = false;
end

%% INITIALIZE
% Check which post processed arrays are available
if isfield(module.Results,'PostProc')
    if isfield(module.Results.PostProc,'PriorPred')
        priorPred_avail = true;
    else
        priorPred_avail = false;
    end
    if isfield(module.Results.PostProc,'PostPred')
        postPred_avail = true;
    else
        postPred_avail = false;
    end
    if isfield(module.Results.PostProc,'PostSample')
        procPostSample_avail = true;
        % determine sample size
        [nIter,nDim,nChains] = size(module.Results.PostProc.PostSample);
    else
        procPostSample_avail = false;
        % determine sample size
        [nIter,nDim,nChains] = size(module.Results.Sample);
    end
    if isfield(module.Results.PostProc,'PriorSample')
        procPriorSample_avail = true;
    else
        procPriorSample_avail = false;
    end
    if isfield(module.Results.PostProc,'PointEstimate')
        pointEstimate_avail = true;
        pointEstimate = module.Results.PostProc.PointEstimate;
        if isfield(module.Results.PostProc,'PostPred')
            if isfield(module.Results.PostProc.PostPred.model,'pointEstimateRun')
                pointEstimatePred_avail = true;
                for ii = 1:module.Internal.nForwardModels
                    pointEstimatePred(ii).Y = module.Results.PostProc.PostPred.model(ii).pointEstimateRun;
                    pointEstimatePred(ii).Type = module.Results.PostProc.PointEstimate.Type;
                end
            else
                pointEstimatePred_avail = false;
            end
        else
            pointEstimatePred_avail = false;
        end
    else
        pointEstimate_avail = false;
    end
else
    %nothing is available
    priorPred_avail = false;
    postPred_avail = false;
    procPostSample_avail = false;
    pointEstimate_avail = false;
    procPriorSample_avail = false;
    
    % determine sample size
    [nIter,nDim,nChains] = size(module.Results.Sample);
end

% get labels for parameters
for ii = 1:length(module.Internal.FullPrior.Marginals)
    currLabel = module.Internal.FullPrior.Marginals(ii).Name;
    % assign to container
    paramLabels{ii} = currLabel;
end

%% Default behavior
% acceptance ratio
Default.plotAcceptance_flag = 'false';
% scatterplot
Default.scatterplotParams = 1:nDim; %all
% predictive distributions
if or(priorPred_avail,postPred_avail)
    %if predictive samples are available plot as well
    Default.plotPredDist_flag = true;
else
    Default.plotPredDist_flag = 'false';
end
% mean convergence
Default.plotMeanConvergence_flag = 'false';
% trace plots
Default.plotTrace_flag = 'false';

%% Check for input arguments
%set optional arguments
if nargin > 1
    % vargin given
    parse_keys = {'acceptance','scatterplot','preddist',...
        'meanconvergence','trace'};
    parse_types = {'p','p','p','p','p'};
    % make NAME lower case
    varargin(1:2:end) = lower(varargin(1:2:end));
    [uq_cline, ~] = uq_simple_parser(varargin, parse_keys, parse_types);
else
    % no varargin, use default options
    uq_cline{1} = Default.plotAcceptance_flag;
    uq_cline{2} = Default.scatterplotParams;
    uq_cline{3} = Default.plotPredDist_flag;
    uq_cline{4} = Default.plotMeanConvergence_flag;
    uq_cline{5} = Default.plotTrace_flag;
end

% 'acceptance' plots the acceptance rate for each chain
if ~strcmp(uq_cline{1}, 'false')
    plotAcceptance_flag = uq_cline{1};
else
    plotAcceptance_flag = false;
end

% 'scatterplot' plots an mDim parameter scatter
if ~strcmp(uq_cline{2}, 'false')
    plotScatterplot_flag = true;
    if isnumeric(uq_cline{2})
        scatterplotParams = uq_cline{2};
    elseif strcmp(uq_cline{2}, 'all')
        scatterplotParams = 1:nDim;
    elseif strcmp(uq_cline{2}, 'none')
        plotScatterplot_flag = false;
    else
        error('Wrong value found in scatterplot name value pair')
    end
else
    plotScatterplot_flag = false;
end

% 'predDist' plots the prior and posterior predictive
% distributions (if available)
if ~strcmp(uq_cline{3}, 'false')
    if ~or(postPred_avail,priorPred_avail)
        % neither prior nor posterior predictive is available
        error('Need to provide prior or posterior predictive model evaluations. See uq_postProcessInversion().')
    end
    if CUSTOM_LIKELI
        error('Predictive distributions are not supported with user-specified likelihood functions.')
    end
    if priorPred_avail
        plotPredDist_flag = true;
    end
    if postPred_avail
        plotPredDist_flag = true;
    end
else
    plotPredDist_flag = false;
end

% 'meanConvergence' plots the convergence of all chains
if ~strcmp(uq_cline{4}, 'false')
    plotMeanConvergence_flag = true;
    if isnumeric(uq_cline{4})
        plotMeanConvergenceIndex = uq_cline{4};
    elseif strcmp(uq_cline{4}, 'all')
        plotMeanConvergenceIndex = 1:nDim;
    elseif strcmp(uq_cline{4}, 'none')
        plotMeanConvergence_flag = false;
    else
        error('Wrong value found in meanConvergence name value pair')
    end
else
    plotMeanConvergence_flag = false;
end

% 'plotTrace' plots the chain trace plots
if ~strcmp(uq_cline{5}, 'false')
    plotTrace_flag = true;
    if isnumeric(uq_cline{5})
        plotTraceIndex = uq_cline{5};
    elseif strcmp(uq_cline{5}, 'all')
        plotTraceIndex = 1:nDim;
    elseif strcmp(uq_cline{5}, 'none')
        plotTrace_flag = false;
    else
        error('Wrong value found in trace name value pair')
    end
else
    plotTrace_flag = false;
end

%% Create the plots

%% Plot the acceptance rate for each chain
if plotAcceptance_flag
    % Retrieve acceptance rates
    Acceptance = module.Results.Acceptance;
    % Open a new figure
    uq_figure
    % Create a scatterplot of the acceptance rates
    uq_plot(1:nChains, Acceptance, 'x', 'MarkerSize', 10)
    ylim([0 1])
    xlabel('Chain')
    ylabel('Acceptance Rate')
    title('Acceptance Rate per Chain')
end

%% Plot the m-dimensional parameter scatterplot
if plotScatterplot_flag
    % Define custom colors
    priorColor = [0.5 0.7 1.0];
    postColor  = [0.0 0.2 0.6];
    
    % Update scatterplot labels
    paramLabelsScatter = {paramLabels{scatterplotParams}};
    
    % If prior sample is available, create the prior plot
    if procPriorSample_avail  
        scatterplotPrior_Sample = reshape(...
            permute(module.Results.PostProc.PriorSample,[2 1 3]),nDim,[])';
        % Get relevant subset of Sample
        PriorSample = scatterplotPrior_Sample(:,scatterplotParams);
        % Open a new figure
        uq_figure
        % Create the scatterplot
        uq_scatterDensity(...
            gca,...
            PriorSample,...
            'Labels', paramLabelsScatter,...
            'Color', priorColor,...
            'Title', 'Prior Sample');
    end
    
    % Check for posterior sample
    if procPostSample_avail  
        scatterplotPost_Sample = reshape(...
            permute(module.Results.PostProc.PostSample,[2 1 3]),nDim,[])';
    else
        scatterplotPost_Sample = reshape(...
            permute(module.Results.Sample,[2 1 3]),nDim,[])';
    end
    % Get relevant subset (requested dimension) of the posterior sample
    PostSample = scatterplotPost_Sample(:,scatterplotParams);
    % Open a new figure
    uq_figure
    % Switch between point estimate and no point estimate
    if pointEstimate_avail
        uq_scatterDensity(...
            gca,...
            PostSample,...
            'Labels', paramLabelsScatter,...
            'Color', postColor,...
            'Title', 'Posterior Sample',...
            'Points', pointEstimate.X(:,scatterplotParams));
    else
        uq_scatterDensity(...
            gca,...
            PostSample,...
            'Labels', paramLabelsScatter,...
            'Color', postColor,...
            'Title', 'Posterior Sample');
    end
end

if plotPredDist_flag
    % Plot samples from the prior and posterior predictive distribution
    % (if available)
    %check which model evaluations are available
    if and(priorPred_avail, postPred_avail)
        %plot prior and posterior predictive (combine both runs)
        model = module.Results.PostProc.PriorPred.model;
        for ii = 1:length(model)
            model(ii).postPredRuns = module.Results.PostProc.PostPred.model(ii).postPredRuns;
        end
        plotType = 'priorPost';
    elseif priorPred_avail
        %plot only prior pred
        model = module.Results.PostProc.PriorPred.model;
        plotType = 'prior';
    elseif postPred_avail
        %plot only posterior pred
        model = module.Results.PostProc.PostPred.model;
        plotType = 'post';
    end
     
    %loop over forward models
    for ii = 1:module.Internal.nForwardModels
        %plot
        if and(size(model(ii).data) == 1,size(model(ii).data(1).y,2) == 1)
            % histogram for scalar model outputs
            if pointEstimatePred_avail
                plotSinglePred(model(ii),plotType,pointEstimatePred(ii))
            else
                plotSinglePred(model(ii),plotType)
            end 
        else
            % line plots for vectorized model outputs
            if pointEstimatePred_avail
                plotSeriesPred(model(ii),plotType,pointEstimatePred(ii))
            else
                plotSeriesPred(model(ii),plotType)
            end
            xlim([0,size(model(ii).data,2)+1])
        end
        %add title to plot
        title(sprintf('Forward model %i',ii),'Fontsize',16)
    end
end

if plotMeanConvergence_flag
    % check for posterior sample
    if procPostSample_avail  
        meanConvergence_Sample = module.Results.PostProc.PostSample;
    else
        meanConvergence_Sample = module.Results.Sample;
    end
    plotIndex = plotMeanConvergenceIndex;
    
    %loop over plot indices
    for ii = plotIndex
        % plot only a certain amount of steps
        Nplotsteps = min(1000,nIter);
        plotSteps = unique(floor(linspace(1,nIter,Nplotsteps)));

        %compute means
        meanVals = zeros(Nplotsteps,nChains);
        for jj = 1:numel(plotSteps)
            currPlotStep = plotSteps(jj);
            % loop over chains
            for kk = 1:nChains
                meanVals(jj,kk) = mean(meanConvergence_Sample(1:currPlotStep,ii,kk)); 
            end
        end
        
        % combine chains
        meanComb = mean(meanVals,2);

        %open figure
        uq_figure
        hold on
        
        uq_plot(plotSteps,meanComb)
        xlabel('Step','Interpreter','latex')
        ylabel(strcat('$E[',paramLabels{ii},']$'),'Interpreter','latex')
    end
end

if plotTrace_flag
    % update traceplot labels 
    paramLabelsScatter = {paramLabels{plotTraceIndex}};
    
    % get relevant sample
    traceSamples = module.Results.Sample(:,plotTraceIndex,:);
    
    % call uq_traceplot function
    uq_traceplot(traceSamples,'labels',paramLabelsScatter);   
end

end

%% ------------------------------------------------------------------------
function plotSinglePred(modelEvals, plotType, varargin)
%PLOTSINGLEPRED creates a plot of single prediction.

% Use histogram for simple data
if nargin > 2
    % With point estimate
    pointEstimatePredFlag = true;
    pointEstimatePred = varargin{1};
else
    % Without point estimate
    pointEstimatePredFlag = false;
end
% Get the data 
y = modelEvals.data.y;

% Use custom colors
priorColor = [0.5, 0.7, 1.0];
postColor  = [0.0, 0.2, 0.6];

%% Create the plot
uq_figure  % Open a new figure
switch plotType
    case 'priorPost'
        % Create histogram plots for both prior & posterior predictive runs
        priorRuns = modelEvals.priorPredRuns;
        postRuns = modelEvals.postPredRuns;
        
        % Plot for prior predictive runs
        priorPlot = uq_histogram(priorRuns, 'Facecolor', priorColor);
        xData = get(priorPlot,'XData');
        
        % Posterior predictive runs
        % Compute the center of the bins of the histogram
        barSpacing = xData(2) - xData(1);
        % If postRuns outside xData, extend xData
        if max(postRuns) > max(xData)
            xData = [xData(end-1) xData(end):barSpacing:max(postRuns)];
        end
        if min(postRuns) < min(xData)
            xData = [min(postRuns):barSpacing:xData(1) xData(2:end)];
        end
        % Plot posterior predictive runs
        postPlot = uq_histogram(postRuns, xData, 'FaceColor', postColor);
        
        % Update the legend
        legendObj = [priorPlot(1) postPlot(1)];
        legendName = {'prior predictive','posterior predictive'};
    
    case 'prior'
        % Create a histogram plot only for prior predictive runs
        priorRuns = modelEvals.priorPredRuns;
        % Plot prior predictive runs
        priorPlot = uq_histogram(priorRuns, 'FaceColor', priorColor);
        % Update the legend
        legendObj = priorPlot(1);
        legendName = {'prior predictive'};
    
    case 'post'
        % Create a histogram plot for posterior predictive runs
        postRuns = modelEvals.postPredRuns;
        % Plot posterior predictive runs
        postPlot = uq_histogram(...
            postRuns,...
            'FaceColor', postColor,...
            'EdgeColor', 'none');
        % Update the legend
        legendObj = postPlot(1);
        legendName = {'posterior predictive'};
        hold on
end

% Add point estimate
if pointEstimatePredFlag
    pointEstimatePlot = uq_plot([pointEstimatePred.Y pointEstimatePred.Y],ylim,'r-');
    % Update legend
    legendObj = [legendObj,pointEstimatePlot(1)];
    legendName{end+1} = lower(pointEstimatePred.Type);
end

% Plot the histogram of the data
dataPlot = uq_plot([y y]',[zeros(size(y)) ones(size(y))]'*0.15*range(ylim),'g');
hold off

% Create the legend
legendObj = [legendObj,dataPlot(1)];
legendName{end+1} = 'data';
uq_legend(legendObj, legendName, 'Location', 'NorthWest')

% Put the axes label
xlabel('$\mathcal{Y}$')

end

%% ------------------------------------------------------------------------
function plotSeriesPred(modelEvals, plotType, varargin)
%PLOTSERIESPRED creates a plot of serial (non-scalar) predictions.

%% Parse inputs
if nargin > 2
    % Point estimate available
    pointEstimatePredFlag = true;
    pointEstimatePred = varargin{1};
else
    % Point estimate not available
    pointEstimatePredFlag = false;
end

% Get the (observed) data
data = modelEvals.data;
xDummy = 1:numel(data);

% Set custom colors
priorColor = [0.5, 0.7, 1.0];
postColor  = [0.0, 0.2, 0.6];

%% Create the plot
uq_figure  % Open a new figure

% Use violin plot for histogram of data series
switch plotType
    case 'priorPost'
        % Both prior and posterior predictive runs
        priorRuns = modelEvals.priorPredRuns;
        postRuns = modelEvals.postPredRuns;
        % Plot prior predictive runs
        priorPlot = uq_violinplot(priorRuns, 'FaceColor', priorColor);
        hold on
        % Plot posterior predictive runs
        postPlot = uq_violinplot(postRuns, 'FaceColor', postColor);
        % Update legend target
        legendObj = [priorPlot(1) postPlot(1)];
        legendName = {'prior predictive', 'posterior predictive'};
    case 'prior'
        % Only prior prior predictive runs
        priorRuns = modelEvals.priorPredRuns;
        % Plot prior predictive runs
        priorPlot = uq_violinplot(priorRuns, 'FaceColor', priorColor);
        % Update legend target
        legendObj = priorPlot(1);
        legendName = {'prior predictive'};
        hold on
    case 'post'
        % Only posterior predictive runs
        postRuns = modelEvals.postPredRuns;
        % Plot posterior predictive runs
        postPlot = uq_violinplot(postRuns, 'FaceColor', postColor);
        % Update legend target
        legendObj = postPlot(1);
        legendName = {'posterior predictive'};
        hold on
end

% Use scatter plot for observed data
for ii = xDummy
    yCurr = data(ii).y;
    for jj = 1:numel(yCurr) 
        dataPlot = scatter(ii, yCurr(jj), 100, 'gx'); 
    end
end

% Update legend target
legendObj = [legendObj dataPlot(1)];
legendName{end+1} = 'data';

% Add point estimate
if pointEstimatePredFlag
    pointEstimatePlot = scatter(xDummy, pointEstimatePred.Y, 100, 'r+');
    % Update legend target
    legendObj = [legendObj,pointEstimatePlot(1)];
    legendName{end+1} = lower(pointEstimatePred.Type);
end
hold off

% Add legend
uq_legend(legendObj,legendName)

% Add axes labels
ylabel('$\mathcal{Y}$', 'Interpreter', 'LaTeX')
xlabel('$\mathrm{Data\,index\,}$(-)', 'Interpreter', 'LaTeX')

% Set xtick labels (maximum 10)
nTicks = 10;
if length(xDummy) > nTicks
    xDummy = ceil(linspace(1, length(xDummy), nTicks));
end
set(gca, 'XTick', xDummy)

% If ticklabel interpreter can be set, update labels to latex
labels = cell(length(xDummy),1);
if isprop(gca,'TickLabelInterpreter')
    for ii = 1:length(xDummy)
        labels{ii} = sprintf('$\\mathrm{y_{%u}}$',xDummy(ii));
    end
else
    for ii = 1:length(xDummy)
        labels{ii} = sprintf('y%d',xDummy(ii));
    end
end
set(gca, 'XTickLabel', labels)

end

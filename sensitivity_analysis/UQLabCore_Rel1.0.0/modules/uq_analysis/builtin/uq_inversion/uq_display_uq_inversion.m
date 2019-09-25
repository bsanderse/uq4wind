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
%                          default : 0
%       'trace'            Plots the trace plot of the marginal sample
%                          points for the specified parameter 
%                          - Integer or 'all'
%                          default : 0
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
        pointParam = module.Results.PostProc.PointEstimate.Parameter;
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
% scatterplot
Default.plotScatterplot_flag = true;
Default.scatterplotParams = 1:nDim; %all
% predictive distributions
if or(priorPred_avail,postPred_avail)
    %if predictive samples are available plot as well
    Default.plotPredDist_flag = true;
else
    Default.plotPredDist_flag = false;
end
% acceptance ratio
Default.plotAcceptance_flag = false;
% mean convergence
Default.plotMeanConvergence_flag = false;
% trace plots
Default.plotTrace_flag = false;

%% Check for input arguments
%set optional arguments
if nargin > 1
    % vargin given
    parse_keys = {'acceptance','scatterplot','predDist',...
        'meanConvergence','trace'};
    parse_types = {'p','p','p','p','p'};
    [uq_cline, ~] = uq_simple_parser(varargin, parse_keys, parse_types);
else
    % no vargin, use default options
    nOpts = 5;
    uq_cline = cell(nOpts,1);
    for ii = 1:nOpts
        uq_cline{ii} = 'false';
    end
end


% 'acceptance' plots the acceptance rate for each chain
if ~strcmp(uq_cline{1}, 'false')
    plotAcceptance_flag = uq_cline{1};
else
    plotAcceptance_flag = Default.plotAcceptance_flag;
end

% 'scatterplot' plots an mDim parameter scatter
if ~strcmp(uq_cline{2}, 'false')
    plotScatterplot_flag = true;
    if isnumeric(uq_cline{2})
        scatterplotParams = uq_cline{2};
    elseif strcmp(uq_cline{2}, 'all')
        scatterplotParams = 1:nDim;
    else
        error('Wrong value found in scatterplot name value pair')
    end
else
    plotScatterplot_flag = Default.plotScatterplot_flag;
    scatterplotParams = Default.scatterplotParams;
end

% 'predDist' plots the prior and posterior predictive
% distributions (if available)
if ~strcmp(uq_cline{3}, 'false')
    if ~or(postPred_avail,priorPred_avail)
        % neither prior nor posterior predictive is available
        error('Need to provide prior or posterior predictive model evaluations. See uq_inversion_postProc().')
    end
    if CUSTOM_LIKELI
        error('Predictive distributions are not supported with user-specified likelihood functions.')
    end
    if priorPred_avail
        plotPredDist_flag = true;
        plotPriorPred_flag = true;
    end
    if postPred_avail
        plotPredDist_flag = true;
        plotPostPred_flag = true;
    end
else
    plotPredDist_flag = Default.plotPredDist_flag;
end

% 'meanConvergence' plots the convergence of all chains
if ~strcmp(uq_cline{4}, 'false')
    plotMeanConvergence_flag = true;
    if isnumeric(uq_cline{4})
        plotMeanConvergenceIndex = uq_cline{4};
    elseif strcmp(uq_cline{4}, 'all')
        plotMeanConvergenceIndex = 1:nDim;
    else
        error('Wrong value found in meanConvergence name value pair')
    end
else
    plotMeanConvergence_flag = Default.plotMeanConvergence_flag;
end

% 'plotTrace' plots the chain trace plots
if ~strcmp(uq_cline{5}, 'false')
    plotTrace_flag = true;
    if isnumeric(uq_cline{5})
        plotTraceIndex = uq_cline{5};
    elseif strcmp(uq_cline{5}, 'all')
        plotTraceIndex = 1:nDim;
    else
        error('Wrong value found in trace name value pair')
    end
else
    plotTrace_flag = Default.plotTrace_flag;
end

%% Plots
fs = 12;

% plot the acceptance rate for each chain
if plotAcceptance_flag
    %retreive acceptance rates
    Acceptance = module.Results.Acceptance;

    %open figure
    figure('pos',[10 10 600 200],'color','white')

    %scatter the Acceptance rate
    scatter(1:nChains,Acceptance)
    ylim([0 1])
    uq_setInterpreters(gca)
    xlabel('Chain', 'Fontsize', fs+2)
    ylabel('Acceptance Rate', 'Fontsize', fs+2)
    title('Acceptance Rate per Chain', 'Fontsize', fs+2)
    grid on
    box on
    set(gca, 'LineWidth', 2, 'Fontsize', fs)
end

% plot the mDim parameter scatter
if plotScatterplot_flag
    % colors
    priorColor = [0.5, 0.7, 1];
    postColor = [0, 0.2, 0.6];
    
    %if prior sample is available create prior plot
    if procPriorSample_avail  
        scatterplotPrior_Sample = reshape(permute(module.Results.PostProc.PriorSample,[2 1 3]),nDim,[])';
        % get relevant subset of Sample
        PriorSample = scatterplotPrior_Sample(:,scatterplotParams);
        % plot
        plotScatter(PriorSample, paramLabels, priorColor, CUSTOM_LIKELI)
        %add title to plot
        axes('Units','Normal');
        h = title('Prior Sample','Interpreter','latex','Fontsize',fs+2);
        set(gca,'visible','off')
        set(h,'visible','on')
    end
    
    %check for posterior sample
    if procPostSample_avail  
        scatterplotPost_Sample = reshape(permute(module.Results.PostProc.PostSample,[2 1 3]),nDim,[])';
    else
        scatterplotPost_Sample = reshape(permute(module.Results.Sample,[2 1 3]),nDim,[])';
    end
    % get relevant subset of Sample
    PostSample = scatterplotPost_Sample(:,scatterplotParams);
    % switch between point estimate and no point estimate
    if pointEstimate_avail
        plotScatter(PostSample, paramLabels, postColor, CUSTOM_LIKELI, pointParam(:,scatterplotParams))
    else
        plotScatter(PostSample, paramLabels, postColor, CUSTOM_LIKELI)
    end
    
    %add title to plot
    axes('Units','Normal');
    h = title('Posterior Sample','Interpreter','latex','Fontsize',fs+2);
    set(gca,'visible','off')
    set(h,'visible','on')
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
        %open figure
        figure('pos',[110 110 600 600],'color','white')
        hold on
        %plot
        if and(size(model(ii).data) == 1,size(model(ii).data(1).y,2) == 1)
            % histogram for scalar model outputs
            plotSinglePred(model(ii),plotType,fs)
        else
            % line plots for vectorized model outputs
            plotSeriesPred(model(ii),plotType,fs)
            xlim([0,size(model(ii).data,2)+1])
        end
        %add title to plot
        uq_setInterpreters(gca)
        title(sprintf('Forward model %i',ii),'Fontsize',fs+2)
        grid on
        box on
        set(gca,'LineWidth', 2, 'Fontsize', fs)
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
        figure('pos',[110 110 600 600],'color','white')
        hold on
        
        uq_setInterpreters(gca)
        xlabel('Step', 'Fontsize', fs)
        ylabel(strcat('$E[',paramLabels{ii},']$'), 'Fontsize',fs)
        plot(plotSteps,meanComb)
        grid on
        box on
        set(gca,'LineWidth',2,'Fontsize',fs)
    end
end

if plotTrace_flag
    % get relevant sample
    chains_Sample = module.Results.Sample;
    % plot the chains specifiec by plotChainsIndex
    nIter = size(chains_Sample,1);
    plotIndex = plotTraceIndex;
    
    %loop over plot indices
    for ii = plotIndex
        %open figure
        figure('pos',[110 110 600 600],'color','white')
        hold on
        for jj = 1:nChains
            plot(1:nIter,chains_Sample(:,ii,jj))   
        end
        uq_setInterpreters(gca)
        xlabel('Step','Fontsize',fs)
        ylabel(paramLabels{ii},'Fontsize',fs)
        grid on
        box on
        set(gca, 'LineWidth', 2, 'Fontsize', fs)
    end          
end

end

function plotSinglePred(modelEvals,plotType,fs)
% histogram for simple data
if isfield(modelEvals,'pointEstimateRun')
    pointEstimate_flag = true;
    pointEstimateRun = modelEvals.pointEstimateRun;
else
    pointEstimate_flag = false;
end
%data is always there
y = modelEvals.data.y;

% colors
priorColor = [0.5, 0.7, 1];
postColor = [0, 0.2, 0.6];
    
switch plotType
    case 'priorPost' %plot both
        priorRuns = modelEvals.priorPredRuns;
        postRuns = modelEvals.postPredRuns;
        %prior runs
        priorPlot = uq_histogram(priorRuns,'Facecolor',priorColor,'edgecolor','none');
        xData = get(priorPlot,'XData');
        %posterior runs
        % if postRuns outside xData, extend xData
        barSpacing = xData(2) - xData(1);
        if max(postRuns) > max(xData)
            xData = [xData(end-1),xData(end):barSpacing:max(postRuns)];
        end
        if min(postRuns) < min(xData)
            xData = [min(postRuns):barSpacing:xData(1),xData(2:end)];
        end
        % plot posterior runs
        postPlot = uq_histogram(postRuns,xData,'Facecolor',postColor,'edgecolor','none');
        %setup legend
        legendObj = [priorPlot(1),postPlot(1)];
        legendName = {'prior predictive','posterior predictive'};
    case 'prior'
        priorRuns = modelEvals.priorPredRuns;
        %prior runs
        priorPlot = uq_histogram(priorRuns,'Facecolor',priorColor,'edgecolor','none');
        %setup legend
        legendObj = priorPlot(1);
        legendName = {'prior predictive'};
    case 'post'
        postRuns = modelEvals.postPredRuns;
        %posterior runs
        postPlot = uq_histogram(postRuns,'Facecolor',postColor,'edgecolor','none');
        %setup legend
        legendObj = postPlot(1);
        legendName = {'posterior predictive'};
end

grid on
box on
uq_setInterpreters(gca)
set(gca,'LineWidth',2,'Fontsize',fs)

%add point estimate
if pointEstimate_flag
    pointEstimatePlot = plot([pointEstimateRun pointEstimateRun],ylim,'r-','linewidth',2);
    %setup legend
    legendObj = [legendObj,pointEstimatePlot(1)];
    legendName{end+1} = 'point estimate';
end

%data
dataPlot = plot([y y]',[zeros(size(y)) ones(size(y))]'*0.15*range(ylim),'g','linewidth',2);
grid on
box on
uq_setInterpreters(gca)
set(gca,'LineWidth',2,'Fontsize',fs)

%setup legend
legendObj = [legendObj,dataPlot(1)];
legendName{end+1} = 'data';

%label
xlabel('$\mathcal{Y}$','Interpreter','Latex','Fontsize',fs+2)
l = legend(legendObj,legendName,'Location','NorthWest');
set(l,'Fontsize',fs,'Interpreter','latex')
end

function plotSeriesPred(modelEvals,plotType,fs)
% violin plots for data series
if isfield(modelEvals,'pointEstimateRun')
    pointEstimate_flag = true;
    pointEstimateRun = modelEvals.pointEstimateRun;
else
    pointEstimate_flag = false;
end
% data is always there
data = modelEvals.data;
xDummy = 1:numel(data);

% colors
priorColor = [0.5, 0.7, 1];
postColor = [0, 0.2, 0.6];

switch plotType
    case 'priorPost' %plot both
        priorRuns = modelEvals.priorPredRuns;
        postRuns = modelEvals.postPredRuns;
        %prior runs
        priorPlot = plotViolins(priorRuns,priorColor);
        %posterior runs
        postPlot = plotViolins(postRuns,postColor);
        %setup legend
        legendObj = [priorPlot(1),postPlot(1)];
        legendName = {'prior predictive','posterior predictive'};
    case 'prior'
        priorRuns = modelEvals.priorPredRuns;
        %prior runs
        priorPlot = plotViolins(priorRuns,priorColor);
        %setup legend
        legendObj = priorPlot(1);
        legendName = {'prior predictive'};
    case 'post'
        postRuns = modelEvals.postPredRuns;
        %posterior runs
        postPlot = plotViolins(postRuns,postColor);
        %setup legend
        legendObj = postPlot(1);
        legendName = {'posterior predictive'};
end

%data
for ii = xDummy
    yCurr = data(ii).y;
    for jj = 1:numel(yCurr) 
        dataPlot = scatter(ii,yCurr(jj),100,'gx'); 
    end
end
%setup legend
legendObj = [legendObj,dataPlot(1)];
legendName{end+1} = 'data';

%add point estimate
if pointEstimate_flag
    pointEstimatePlot = scatter(xDummy,pointEstimateRun,100,'r+');
    %setup legend
    legendObj = [legendObj,pointEstimatePlot(1)];
    legendName{end+1} = 'point estimate';
end

% label
ylabel('$\mathcal{Y}$','Interpreter','Latex')
xlabel('$\mathrm{Data\,index\,}$(-)', 'Interpreter', 'latex')
l = legend(legendObj,legendName,'Location','best');
set(l,'Fontsize',fs,'Interpreter','latex')
% set xtick labels (maximum 10)
nTicks = 10;
if length(xDummy) > nTicks
    xDummy = ceil(linspace(1,length(xDummy),nTicks));
end
set(gca,'XTick',xDummy)
% if ticklabel interpreter can be set update labels to latex
if isfield(get(gca),'TickLabelInterpreter')
    for ii = 1:length(xDummy)
        label{ii} = sprintf('$\\mathrm{y_{%u}}$',xDummy(ii));
    end
    set(gca,'XTickLabel',label)
end

grid on
box on
uq_setInterpreters(gca)
set(gca, 'LineWidth', 2, 'Fontsize', fs)
end

function plotScatter(Sample, paramLabels, color, CUSTOM_LIKELI, pointParam)
% initialize
nDim = size(Sample,2);
fs = 12;

% open figure
scatterFig = figure('pos',[110 110 600 600],'color','white');

% get smalles and largest samples in each dimension to set axis properly
if nargin > 4
    minSample = min([Sample;pointParam]);
    maxSample = max([Sample;pointParam]);
else
    minSample = min(Sample);
    maxSample = max(Sample);
end

% colormap for plotHistogram2
% map between color and yellow
map = [linspace(color(1),1,64)',linspace(color(2),1,64)',linspace(color(3),0,64)'];

%loop over variable combinations to plot
for ii = 1:nDim
    for jj = 1:nDim
        %switch between histogram and scatter
        if ii == jj %histogram
            scatterSub(ii) = subplotInv(nDim,nDim,(ii-1)*nDim+jj);
            % plot histogram with scott rule binning
            plotHistogram(Sample(:,ii),color)
            if or(~(jj == 1),ii == 1)
                set(gca,'YTickLabel',[])
            end
            if ~(ii == nDim)
                set(gca,'XTickLabel',[])
            end
            %limits
            xlim([minSample(ii),maxSample(ii)])
            %fix aspect ratio
            pbaspect([1 1 1])

            %display point estimate also
            if nargin > 4
                hold on
                pointEstLine(ii) = plot([pointParam(:,ii),pointParam(:,ii)],get(gca,'ylim'),'r','linewidth',2);
                if isfield(scatterFig,'SizeChangedFcn')
                    % add callback to figure
                    set(scatterFig,'SizeChangedFcn',{@resizeui, pointEstLine, scatterSub});
                end
            end
            
            % plot setup
            grid off
            box on
            uq_setInterpreters(gca)
            set(gca,'LineWidth', 1, 'FontSize', fs-4)
            if isfield(scatterFig,'SizeChangedFcn')
                % trigger SizeChangedFcn callback
                notify(gcf,'SizeChanged')
            end
        elseif jj < ii
            subplotInv(nDim,nDim,(ii-1)*nDim+jj);
            % plot
            plotHistogram2(Sample(:,jj),Sample(:,ii),map)
            if ~(jj == 1)
                set(gca,'YTickLabel',[])
            end
            if ~(ii == nDim)
                set(gca,'XTickLabel',[])
            end
            %limits
            ylim([minSample(ii),maxSample(ii)])
            xlim([minSample(jj),maxSample(jj)])
            %fix aspect ratio
            pbaspect([1 1 1])

            %display point estimate also
            if nargin > 4
                hold on
                plot(pointParam(:,jj),pointParam(:,ii),'w+','MarkerSize',8,'LineWidth',2,'LineStyle', 'none')
                plot(pointParam(:,jj),pointParam(:,ii),'r+','MarkerSize',8,'LineWidth',1,'LineStyle', 'none')
            end
            
            % plot setup
            grid off
            box on
            uq_setInterpreters(gca)
            set(gca,'LineWidth',1,'FontSize',fs-4)
        end

        %add axis labels
        if jj == 1
            ylabel(paramLabels{ii},'Interpreter','Latex','Fontsize',fs+2)
            if ii > 1 % not for first subplot
                % number of ticks and format
                NumTicks = 2;
                if which('ytickformat')
                    ytickformat('%.2f')
                end
                if which('ytickangle')
                    ytickangle(45)
                end
                L = get(gca,'YLim'); L(1) = L(1) + diff(L)*0.05; L(2) = L(2) - diff(L)*0.05;
                set(gca,'YTick',linspace(L(1),L(2),NumTicks))
            end
        end
        if ii == nDim
            xlabel(paramLabels{jj}, 'Interpreter','Latex','Fontsize',fs+2)
            % number of ticks and format
            NumTicks = 2;
            if which('xtickformat')
                xtickformat('%.2f')
            end
            if which('xtickangle')
                xtickangle(45)
            end
            L = get(gca,'XLim'); L(1) = L(1) + diff(L)*0.05; L(2) = L(2) - diff(L)*0.05;
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        end
    end
end
end

function plotHistogram2(histData1, histData2,map)
    min1 = min(histData1); max1 = max(histData1); 
    min2 = min(histData2); max2 = max(histData2); 
    % scott rule
    binWidth1 = 3.5*std(histData1(:))*numel(histData1)^(-1/3); 
    nBins1 = ceil(range(histData1)/binWidth1);
    binWidth2 = 3.5*std(histData2(:))*numel(histData2)^(-1/3); 
    nBins2 = ceil(range(histData2)/binWidth2);
    % build array for plotting
    C = zeros(nBins1+1,nBins2+1);
    for ii = 1:length(histData1)
        index1 = ceil((histData1(ii)-min1)/binWidth1)+1;
        index2 = ceil((histData2(ii)-min2)/binWidth2)+1;
        C(index1,index2) = C(index1,index2) + 1;
    end
    % scale C
    nColor = size(map,1);
    C = ceil(C.*(nColor/max(max(C))));
    % plot
    im = imagesc([min1 max1],[min2 max2],C');
    colormap([[1 1 1];map]);
    set(gca,'YDir','normal')
    set(im,'CDataMapping','direct')
end

function plotHistogram(histData,color)
    % scott rule
    binWidth = 3.5*std(histData(:))*numel(histData)^(-1/3); 
    nBins = ceil(range(histData)/binWidth);
    hist(histData,nBins)
    h = findobj(gca,'Type','patch');
    set(h,'LineStyle','none','FaceColor',color,'FaceAlpha',1)
end

function handle = plotViolins(data,color)
for ii = 1:size(data,2)
    currRuns = data(:,ii);
    [f,xi] = ksdensity(currRuns);
    %scale and shift
    patchF = [f -fliplr(f)];
    patchF = patchF/(2*max(patchF));
    patchF = patchF + ii;
    
    %plot
    handle = patch(patchF,[xi fliplr(xi)],color,'EdgeColor','none');
end
end

function handle = subplotInv(NRow,NCol,index)
% custom subplot with reduced margin
xMarginLab = 0.15;
yMarginLab = 0.15;
xMarginTot = 0.05;
yMarginTot = 0.05;
xMarginLoc = xMarginTot/(NCol);
yMarginLoc = yMarginTot/(NCol);
xLength = (1-xMarginLab*2-xMarginTot)/NCol;
yLength = (1-yMarginLab*2-yMarginTot)/NRow;

% x/y index
[xIndex,yIndex] = ind2sub([NRow,NCol],index);

% coordinates
xCoord = xMarginLab + xMarginLoc*(xIndex-1) + xLength*(xIndex-1);
yCoord = 1 - (yMarginLab + yMarginLoc*yIndex + yLength*(yIndex));

% create axes
handle = axes('Position',[xCoord yCoord xLength yLength]);
end

function resizeui(hObject,event,LineContainer, AxesContainer)
% loop over lines and set y limit dynamically if the figure is resized
for ii = 1:length(LineContainer)
    yLimCurr = AxesContainer(ii).YLim;
    set(LineContainer(ii),'YData',yLimCurr);
end
end

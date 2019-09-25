function uq_PCK_display(PCKRGModel, outArray, varargin)
%Plots the mean and variance of a PC-Kriging predictor. Only works for 1- and 
% 2-D inputs.

% number of discretization points for the figures
N1d = 500;
N2d = 80;

% check whether the PC-Kriging model has been calculated
if ~PCKRGModel.Internal.Runtime.isCalculated
    fprintf('PC-Kriging object %s is not yet initialized!\nGiven Configuration Options:', PCKRGModel.Name);
    PCKRGModel.Options
    return;
end


%% Consistency checks and command line parsing
if ~exist('outArray', 'var')
    outArray = 1;
    if length(PCKRGModel.PCK) > 1;
        warning('The selected PC-Kriging metamodel has more than 1 output. Only the 1st output will be printed');
        fprintf('You can specify the outputs you want to be displayed with the syntax:\n')
        fprintf('uq_print(PCKRGModule, OUTARRAY)\nwhere OUTARRAY is the index of desired outputs, e.g. 1:3 for the first three\n\n')
    end
end
if max(outArray) > length(PCKRGModel.PCK)
    error('Requested output range is too large') ;
end


%% parsing the residual command line
% initialization
if nargin > 2
    parse_keys = {'R','nolegend'};
    parse_types = {'f','f'};
    [uq_cline, varargin] = uq_simple_parser(varargin, parse_keys, parse_types);
    
    % 'R' option additionally prints R matrix
    R_flag = strcmp(uq_cline{1}, 'true');
    noLegend_flag = strcmp(uq_cline{2}, 'true');
        
    flagWasSet = R_flag;
else
    flagWasSet = false;
    noLegend_flag = false;
end



%% Produce non-default plots if any flag  was set
if flagWasSet
    if R_flag
        for ii = 1:length(outArray)
            current_output = outArray(ii);
            uq_figure('name',sprintf('Output #%i', current_output), 'Position', [50 50 500 400])
            imagesc(PCKRGModel.Internal.Kriging(current_output).Internal.Kriging.GP.R );
            xlabel('Column index','FontSize', 14,'Interpreter','latex');
            ylabel('Row index','FontSize', 14,'Interpreter','latex');
            title('R matrix values','Interpreter','latex')
        end
    end
    return
end


%% Produce plot
% Get input dimension
M = PCKRGModel.Internal.Runtime.M ;


switch M
    case 1
        %get experimental design
        X = PCKRGModel.ExpDesign.X;
        Y = PCKRGModel.ExpDesign.Y;
        
        %define prediction points
        Xmin = min(X); Xmax = max(X);
        Xval = linspace(Xmin, Xmax, N1d)';
        Uval = uq_GeneralIsopTransform(Xval, PCKRGModel.Internal.Input.Marginals, PCKRGModel.Internal.Input.Copula,...
            PCKRGModel.Internal.AuxSpace.Marginals, PCKRGModel.Internal.AuxSpace.Copula);
        
        %evaluate prediction points
        [Ymu_KRG, Ysigma_KRG]= uq_evalModel(PCKRGModel,Xval);
        confInterval = 0.05;
        Conf = norminv(1 - 0.5*confInterval, 0, 1 )* sqrt(Ysigma_KRG) ;
        
        for ii = 1:length(outArray)
            %determine the response of the trend alone (i.e. the PCE model
            %multiplied with coefficients and finally summed)
            t_ii = PCKRGModel.Internal.Kriging(ii).Internal.Kriging.Trend.Handle;
            Ytrend = sum(ones(N1d,1)*PCKRGModel.PCK.beta' .* t_ii(Uval), 2);
           
            % Draw the figure
            legendHandles = [];
            legendTxt = {};
            current_output = outArray(ii);
            uq_figure('name', sprintf('Output #%i', current_output),...
                'Position', [50 50 500 400])
            h = uq_plotConfidence(Xval, Ymu_KRG(:,current_output),...
                Conf(:,current_output));
            hold on
            h2 = plot(Xval, Ytrend,...
                'k--', 'LineWidth', 2, 'Color', [1 0.5 0]);
            legendHandles = [legendHandles;  h(:); h2(:)];
            legendTxt = [legendTxt 'PC-Kriging approximation'];
            legendTxt = [legendTxt '$95$\% confidence interval'];
            legendTxt = [legendTxt 'PC trend in PC-Kriging'];
            h = uq_plot(X,Y(:,current_output),...
                'ko', 'MarkerFaceColor','k');
            legendHandles = [legendHandles; h(:)];
            legendTxt = [legendTxt 'Observations'];
            hold off
            if ~noLegend_flag
               leg = legend(legendHandles, legendTxt,...
                   'Location', 'Best', 'Interpreter', 'latex');
            end
            uq_setInterpreters(gca)
            xlabel('$\mathrm{X}$')
            ylabel('$\mathrm{Y}$')
            set(gca, 'FontSize', 14)
        end
        
    case 2
        % Get the experimental design
        X = PCKRGModel.ExpDesign.X;
        X1min = min(X(:,1)); X1max = max(X(:,1));
        X2min = min(X(:,2)); X2max = max(X(:,2));
        
        % Compute a regular mesh
        [X1val, X2val] = meshgrid(linspace(X1min, X1max, N2d), ...
            linspace(X2min, X2max, N2d));
        X1val_v = reshape(X1val,[],1);
        X2val_v = reshape(X2val,[],1);
        
        % Evaluate the PC-Kriging model on the mesh points
        Xval = [X1val_v, X2val_v] ;
        [Ymu_KRG_v, Ysigma_KRG_v] = uq_evalModel(PCKRGModel,Xval);
        
        for ii = 1:length(outArray)
            % Draw the figure
            current_output = outArray(ii);
            uq_figure('name',sprintf('Output #%i', current_output),...
                'Position', [50 50 500 400]);
            
            Ymu_KRG = reshape(Ymu_KRG_v(:,outArray(ii)), size(X1val));
            Ysigma_KRG = reshape(Ysigma_KRG_v(:,outArray(ii)), size(X1val));
            caxisMin = min([Ymu_KRG_v(:);Ysigma_KRG_v(:)]);
            caxisMax = max([Ymu_KRG_v(:);Ysigma_KRG_v(:)]);
            
            % Subplot for the prediction mean value
            subplot(1, 2, 1)
            h = pcolor(X1val, X2val, Ymu_KRG);
            set(h, 'EdgeColor', 'none')
            hold on
            uq_plot(X(:,1), X(:,2),...
                'ko', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
            hold off
            axis([X1min X1max X2min X2max])
            
            uq_setInterpreters(gca)
            xlabel('$\mathrm{X_1}$')
            ylabel('$\mathrm{X_2}$')
            title('$\mathrm{\widehat{\mu}_Y}$')
            caxis([caxisMin, caxisMax])
            set(gca, 'Box' , 'on', 'Layer', 'top', 'FontSize', 14)
            set(colorbar, 'Visible', 'off')
            
            % Subplot for the prediction variance
            subplot(1, 2, 2)
            h = pcolor(X1val, X2val, abs(Ysigma_KRG));
            set(h, 'EdgeColor', 'none');
            hold on
            uq_plot(X(:,1),X(:,2), 'ko', 'MarkerFaceColor','r',...
                'MarkerEdgeColor','none');
            hold off
            axis([X1min X1max X2min X2max])
            
            uq_setInterpreters(gca)
            xlabel('$\mathrm{X_1}$');
            ylabel('$\mathrm{X_2}$');
            title('$\mathrm{\widehat{\sigma}_Y^2}$')
            set(gca, 'Box' , 'on', 'Layer', 'top', 'FontSize', 14)
            caxis([caxisMin, caxisMax])
            colorbar
            
        end
    otherwise
        error('Only 1 and 2 dimensional X''s are supported!')
end
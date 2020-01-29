function uq_SVR_display(SVRModel, outArray, varargin)
%Plots the mean and variance of a Kriging predictor. Only works for 1- and
% 2-D inputs.

%% Internal parameters
% BoundColor = [128,128,128]/255; %gray
BoundColor = [176 196 222]/255; %lightskyblue
% Granularity of the grid in one- and two-D cases in each directioni
N1d = 100;
N2d = 50;

% Check that the tmodel has been computed with success
if ~SVRModel.Internal.Runtime.isCalculated
    fprintf('SVR object %s is not yet calculated!\nGiven Configuration Options:', SVRModel.Name);
    SVRModel.Options
    return;
end

%% Consistency checks
if ~exist('outArray', 'var')
    outArray = 1;
    if length(SVRModel.SVR) > 1
        warning('The selected SVR metamodel has more than 1 output. Only the 1st output will be printed');
        fprintf('You can specify the outputs you want to be displayed with the syntax:\n')
        fprintf('uq_print(SVRModule, OUTARRAY)\nwhere OUTARRAY is the index of desired outputs, e.g. 1:3 for the first three\n\n')
    end
end
if max(outArray) > length(SVRModel.SVR)
    error('Requested output range is too large') ;
end

%% Produce plot
% Get input dimension
M = SVRModel.Internal.Runtime.M ;

switch M
    case 1
        % plot for one-dimensional functions
        for ii = 1:length(outArray)
            % Get relevant output
            current_output = outArray(ii);
            % Get coefficients of the SVR expansion
            beta = SVRModel.SVR(current_output).Coefficients.beta;
            % Get indices of support vectors
            Isv = SVRModel.SVR(current_output).Coefficients.SVidx;
            if SVRModel.Internal.SVR(current_output).OutputScaling == 1
                epsilon = SVRModel.SVR(current_output).Hyperparameters.epsilon ...
                    * SVRModel.Internal.Runtime.stdY(:,current_output) + ...
                    SVRModel.Internal.Runtime.muY(:,current_output);
            else
                epsilon = SVRModel.SVR(current_output).Hyperparameters.epsilon;
            end
            % Get training points
            Xtrain = SVRModel.ExpDesign.X;
            Ytrain = SVRModel.ExpDesign.Y(:,current_output);
            % Define bounds of the plot area
            Xmin = min(Xtrain); Xmax = max(Xtrain);
            % Sample points for plot
            Xval = linspace(Xmin, Xmax, N1d)';
            % Evaluate the SVR model and bounds (w.r.t. to epsilon)
            Ysvr= uq_evalModel(SVRModel,Xval);
            ysvr_upper = Ysvr + epsilon;
            Ysvr_lower = Ysvr - epsilon;
            
            %% Create plot
            legendHandles = [];
            legendTxt = {};
            current_output = outArray(ii);
            uq_figure('name', sprintf('Output #%i', current_output));
            % Plot epsilon insensitive tube
            h = uq_plotConfidence(Xval, Ysvr, epsilon);
            hold on
            % NOTE: 'hold on' command in R2014a does not cycle through
            % the color order, so it must be set manually.
            ax = gca;
            colorOrder = get(ax,'ColorOrder');
            legendHandles = [legendHandles; h(:)];
            legendTxt = [legendTxt 'SVR approximation'];
            legendTxt = [legendTxt 'eps-insensitive tube'];
            % Plot non support vector training points
            Obs = 1:length(Ytrain);
            nIsv = ~ismember(Obs,Isv);
            h1 = uq_plot(Xtrain(nIsv,:), Ytrain(nIsv,current_output), 'ko');
            % Plot support vectors
            h2 = uq_plot(...
                Xtrain(Isv,:), Ytrain(Isv,current_output), 'd',...
                'MarkerFaceColor', colorOrder(2,:),...
                'Color', 'k');
            if isempty(h1)
                legendHandles = [legendHandles; h2(:)];
                legendTxt = [legendTxt 'Support vectors'];
            elseif isempty(h2)
                legendHandles = [legendHandles; h1(:)];
                legendTxt = [legendTxt 'Observations (non SVs)'];
            else
                legendHandles = [legendHandles  ; h1(:); h2(:)] ;
                legendTxt = [legendTxt, 'Observations (non SVs)',...
                    'Support vectors'];
            end
            hold off
            % Put legend
            uq_legend(legendHandles, legendTxt)
            % Set axes labels
            xlabel('$\mathrm{X}$')
            ylabel('$\mathrm{Y}$')

        end
    case 2
        for ii = 1:length(outArray)
            % Get relevant output
            current_output = outArray(ii);
            % Get training points
            Xtrain = SVRModel.ExpDesign.X;
            % Get indices of support vectors
            Isv = SVRModel.SVR(current_output).Coefficients.SVidx;
            % Get coefficients of the expansion
            beta = SVRModel.SVR(current_output).Coefficients.beta;
            % Define bounds of the plot area
            X1min = min(Xtrain(:,1)); X1max = max(Xtrain(:,1));
            X2min = min(Xtrain(:,2)); X2max = max(Xtrain(:,2));
            % Create grid for plot
            [X1val,X2val] = meshgrid(linspace(X1min, X1max, N2d), ...
                linspace(X2min, X2max, N2d));
            X1val_v = reshape(X1val,[],1);
            X2val_v = reshape(X2val,[],1);
            
            Xval = [X1val_v, X2val_v] ;
            % Eval SVR model
            Ysvr = uq_evalModel(SVRModel,Xval);
            Ysvr_grid = reshape(Ysvr, size(X1val));
            % Create plot
            uq_figure('name', sprintf('Output #%i', current_output));
            % Plot SVR prediction
            h = pcolor(X1val, X2val, Ysvr_grid)
            hold on
            shading interp
            uq_formatDefaultAxes(gca);
            % Plot training points
            h1 = uq_plot(Xtrain(~Isv,1), Xtrain(~Isv,2),...
                'ko', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
            h2 = uq_plot(Xtrain(Isv,1), Xtrain(Isv,2),...
                'kd', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none');
            cbar = colorbar('Location', 'eastoutside');
            hold off
            % Further formatting of the plot
            axis([X1min X1max X2min X2max])
            % Set labels
            xlabel('$\mathrm{X_1}$')
            ylabel('$\mathrm{X_2}$')
        end
    otherwise
        error('Only 1 and 2 dimensional X''s are supported!')
end

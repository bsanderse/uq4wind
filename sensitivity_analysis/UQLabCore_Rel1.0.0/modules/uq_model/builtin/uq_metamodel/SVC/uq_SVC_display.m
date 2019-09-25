function uq_SVC_display(SVCModel, outArray, varargin)
% UQ_SVC_DISPLAY(SVCMODEL,OUTARRAY,VARARGIN): plot the classes and
% margins specified by the SVC model
% SVCMODEL. Only works for 2-D inputs.
%
% See also: UQ_SVR_DISPLAY, UQ_DISPLAY_UQ_METAMODEL

%% Some internal parameters
% Granularity of the grid in each direction for the 2D plot
granul = 100 ;

%%
% Check that the model has been computed with success
if ~SVCModel.Internal.Runtime.isCalculated
    fprintf('SVC object %s is not yet calculated!\nGiven Configuration Options:', SVCModel.Name);
    SVCModel.Options
    return;
end


%% parse varargin options
% initialization
legend_flag = false;

if nargin > 2
    parse_keys = {'legend'};
    parse_types = {'f', 'f'};
    [uq_cline, varargin] = uq_simple_parser(varargin, parse_keys, parse_types);
    % 'coefficients' option additionally prints the coefficients
    if strcmp(uq_cline{1}, 'true')
        legend_flag = true;
    end
end


%% Consistency checks
if ~exist('outArray', 'var')
    outArray = 1;
    if length(SVCModel.SVC) > 1
        warning('The selected SVC metamodel has more than 1 output. Only the 1st output will be printed');
        fprintf('You can specify the outputs you want to be displayed with the syntax:\n')
        fprintf('uq_print(SVCModule, OUTARRAY)\nwhere OUTARRAY is the index of desired outputs, e.g. 1:3 for the first three\n\n')
    end
end
if max(outArray) > length(SVCModel.SVC)
    error('Requested output range is too large') ;
end

%% Produce plot
% Get input dimension
M = SVCModel.Internal.Runtime.M;

if M == 2
    % Get experimental design input
    Xtrain = SVCModel.ExpDesign.X;
    % Define bounds of the plot area
    Xmin = min(Xtrain);
    Xmax = max(Xtrain);
    minX1 = floor(Xmin(1));
    maxX1 = ceil(Xmax(1));
    minX2 = floor(Xmin(2));
    maxX2 = ceil(Xmax(2));
    % Define grid on the plot area
    [X1grid,X2grid] = meshgrid(linspace(minX1,maxX1,granul),...
        linspace(minX2,maxX2,granul));
    mx = [X1grid(:) X2grid(:)];
    % Evaluate the SVC model
    [Yclass_all,Yval_all] = uq_evalModel(SVCModel,mx);
    for ii = 1:length(outArray)
        % Get desired output
        current_output = outArray(ii);
        Isv = SVCModel.SVC(current_output).Coefficients.SVidx;
        % Get experimental design output
        Ytrain = SVCModel.ExpDesign.Y(:,current_output);
        % Get support vectors
        Xsv = Xtrain(Isv,:);
        % Current output prediction
        Yval = Yval_all(:,current_output) ;
        % Create grid
        Yval_grid = reshape(Yval', granul, granul);
        
        %% Start plotting
        uq_figure('Position',[50 50 500 400])
        % Plot the training sets - except the support vectors that will be
        % added later
        TrnMarker_1 = 'o';
        TrnMarker_2 = 'square';
        TrnMarkersize = 6; 
        TrnLineWidth = 1;
        % Get support vector (is support vector?)
        isSV = zeros(length(Ytrain),1);
        isSV(Isv) = 1;
        % Non-support vectors of the first class
        isSV_1 = ~isSV .* Ytrain > 0;
        % Non-support vectors of the second class
        isSV_2 = ~isSV .* Ytrain < 0;
        h1 = plot(Xtrain(isSV_1,1), Xtrain(isSV_1,2), TrnMarker_1);
        hold on
        set(h1, 'Markersize', TrnMarkersize, 'LineWidth', TrnLineWidth,...
            'Color', 'blue', 'MarkerFaceColor', 'blue');
        h2 = plot(Xtrain(isSV_2,1), Xtrain(isSV_2,2), TrnMarker_2);
        set(h2, 'Markersize', TrnMarkersize, 'LineWidth', TrnLineWidth,...
            'Color', 'red', 'MarkerFaceColor', 'red');
        % Mark the support vectors with green diamonds.
        % (The size of each is not anymore computed relatively
        % w.r.t. the SVC coefficient. - Release 1.1.0)
        h = plot(Xsv(:,1), Xsv(:,2), 'kd', 'MarkerSize', 7,...
            'MarkerFaceColor', 'g');
        % Plot separating line : isoline x \in \Xx : M_{svc}(x) == 0
        [c,hs] = contour(X1grid, X2grid, Yval_grid, [0 0],...
            'k', 'LineWidth', 4);
        % Plot lower margin : isoline x \in \Xx : M_{svc}(x) == -1
        [c,hr] = contour(X1grid, X2grid, Yval_grid, [-1 -1],...
            '--', 'Color', 'k', 'LineWidth', 2);
        % Plot upper margin : isoline x \in \Xx : M_{svc}(x) == 1
        [c,hb] = contour(X1grid, X2grid, Yval_grid, [1 1],...
            '--', 'Color', 'k', 'LineWidth', 2);
        grid on
        % Draw a legend
        if legend_flag
            set(gcf,'Position', [50 50 700 400])
            lg = legend([h(1) h1 h2 hs hr],...
                {'Support vectors' 'Train. set - Cat.1'...
                'Train. set - Cat.2' 'Classifier' 'Margins'});
            set(lg, 'Location', 'eastoutside', 'Interpreter', 'latex')
            legend boxoff
        end
        % Format the resulting plot
        axis([ minX1 maxX1 minX2 maxX2])
        uq_setInterpreters(gca)
        xlabel('$\mathrm{X_1}$')
        ylabel('$\mathrm{X_2}$')
        set(findall(gcf, '-property', 'FontSize'), 'FontSize', 14)
        set(gca, 'LineWidth', 2, 'Box', 'on')
        hold off
    end
else
    error('Only two-dimensional X''s are supported!')
end

end

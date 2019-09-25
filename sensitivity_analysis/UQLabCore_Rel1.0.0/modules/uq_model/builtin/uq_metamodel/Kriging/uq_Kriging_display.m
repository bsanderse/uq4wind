function uq_Kriging_display(KRGModel, outArray, varargin)
% UQ_KRIGING_DISPLAY(KRGMODEL,OUTARRAY,VARARGIN): plot the mean and
% variance of the Kriging predictor specified by the Kriging model
% KRGMODEL. Only works for 1- and 2-D inputs. 
%
% See also: UQ_PCE_DISPLAY, UQ_DISPLAY_UQ_METAMODEL

%% Parameters
N1d = 500;
N2d = 80;

%% Consistency checks and command line parsing
if ~KRGModel.Internal.Runtime.isCalculated
    fprintf('Kriging object %s is not yet initialized!\nGiven Configuration Options:', KRGModel.Name);
    KRGModel.Options
    return;
end

if ~exist('outArray', 'var')
    outArray = 1;
    if length(KRGModel.Kriging) > 1;
        warning('The selected Kriging metamodel has more than 1 output. Only the 1st output will be printed');
        fprintf('You can specify the outputs you want to be displayed with the syntax:\n')
        fprintf('uq_print(KRGModule, OUTARRAY)\nwhere OUTARRAY is the index of desired outputs, e.g. 1:3 for the first three\n\n')
    end
end
if max(outArray) > length(KRGModel.Kriging)
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
            imagesc(KRGModel.Internal.Kriging(current_output).GP.R );
            xlabel('Column index','FontSize', 14,'Interpreter','latex');
            ylabel('Row index','FontSize', 14,'Interpreter','latex');
            title('R matrix values','Interpreter','latex')
        end
    end
    return
end


%% Produce plot
% Get input dimension
M = KRGModel.Internal.Runtime.M ;
switch M
    case 1
        X = KRGModel.ExpDesign.X;
        Y = KRGModel.ExpDesign.Y;
        Xmin = min(X); Xmax = max(X);
        Xval = linspace(Xmin, Xmax, N1d)';
        [Ymu_KRG, Ysigma_KRG]= uq_evalModel(KRGModel,Xval);
        confInterval = 0.05;
        Conf = norminv(1 - 0.5*confInterval, 0, 1 )* sqrt(Ysigma_KRG);
        
        for ii = 1:length(outArray)
            legendHandles = [];
            legendTxt = {};
            current_output = outArray(ii);
            uq_figure('name',sprintf('Output #%i', current_output),...
                'Position', [50 50 500 400])
            h = uq_plotConfidence(Xval, Ymu_KRG(:,current_output),...
                Conf(:,current_output));
            legendHandles = [legendHandles;  h(:)];
            legendTxt = [legendTxt 'Kriging approximation'];
            legendTxt = [legendTxt '$95\%$ confidence interval'];
            hold on
            h = uq_plot(X,Y(:,current_output),...
                'ko', 'MarkerFaceColor', 'k');
            legendHandles = [legendHandles; h(:)];
            legendTxt = [legendTxt 'Observations'];
            hold off
            
            if ~noLegend_flag
                leg = legend(legendHandles, legendTxt);
                set(leg, 'Location', 'best', 'FontSize', 14,...
                    'Interpreter','latex')
            end
            
            uq_setInterpreters(gca)
            xlabel('$\mathrm{X}$')
            ylabel('$\mathrm{Y}$')
            set(gca, 'FontSize', 14)
        end
    case 2
        for ii = 1:length(outArray)
            current_output = outArray(ii);
            
            X = KRGModel.ExpDesign.X;
            
            X1min = min(X(:,1)); X1max = max(X(:,1));
            X2min = min(X(:,2)); X2max = max(X(:,2));
            
            [X1val, X2val] = meshgrid(linspace(X1min, X1max, N2d), ...
                linspace(X2min, X2max, N2d));
            X1val_v = reshape(X1val,[],1);
            X2val_v = reshape(X2val,[],1);
            
            Xval = [X1val_v, X2val_v] ;
            [Ymu_KRG_v, Ysigma_KRG_v] = uq_evalModel(KRGModel,Xval);
            Ymu_KRG = reshape(Ymu_KRG_v, size(X1val));
            Ysigma_KRG = reshape(Ysigma_KRG_v, size(X1val));
            
            % Kriging predictor mean figure
            uq_figure('name',...
                sprintf('Output #%i (mean)', current_output),...
                'Position', [50 50 500 400])

            h = pcolor(X1val, X2val, Ymu_KRG);
            set(h, 'EdgeColor', 'none')
            shading interp
            hold on
            plot(X(:,1), X(:,2),...
                'ko', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none')
            hold off
            axis([X1min X1max X2min X2max])

            cbar = colorbar('Location', 'eastoutside');
            uq_setInterpreters(cbar)

            uq_setInterpreters(gca)
            xlabel('$\mathrm{X_1}$');
            ylabel('$\mathrm{X_2}$');
            title('$\mathrm{\mu_{\widehat{Y}}}$')
            set(gca, 'fontsize', 14);
            
            % Kriging predictor variance figure
            uq_figure('name',...
                sprintf('Output #%i (variance)', current_output),...
                'Position', [50 50 500 400]);
            h = pcolor(X1val, X2val, abs(Ysigma_KRG));
            set(h, 'EdgeColor', 'none')
            shading interp
            hold on
            plot(X(:,1),X(:,2),...
                'ko', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'none')
            hold off
            axis([X1min X1max X2min X2max])

            cbar = colorbar('Location','eastoutside');
            uq_setInterpreters(cbar)

            uq_setInterpreters(gca)
            xlabel('$\mathrm{X_1}$');
            ylabel('$\mathrm{X_2}$');
            title('$\mathrm{\sigma_{\widehat{Y}}^2}$')
            set(gca, 'Box' , 'on', 'Layer', 'top', 'FontSize', 14)
        end

    otherwise
        error('Only 1 and 2 dimensional X''s are supported!')
end

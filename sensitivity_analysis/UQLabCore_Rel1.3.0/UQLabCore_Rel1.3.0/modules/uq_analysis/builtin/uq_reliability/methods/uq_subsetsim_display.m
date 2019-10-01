function varargout = uq_subsetsim_display( module, idx, varargin )
% UQ_SUBSETSIM_DISPLAY visualizes the analysis and results of subset
% simulation
% 
% See also: UQ_DISPLAY_UQ_RELIABILITY

%check whether the model evaluations were stored
if module.Internal.SaveEvaluations
    
    %for each index
    for oo = idx
        
        %display the samples of each subset in 1- and 2-dimensional cases
        switch length(module.Internal.Input.Marginals)
            % Scatter plots of the subset sample for 2-dimensional problems
            case 2
                uq_figure
                hold on
                ax = gca;
                uq_formatDefaultAxes(ax)
                % NOTE: 'hold on' command in R2014a does not cycle through
                % the color order so it must be set manually. 
                colorOrder = get(ax,'ColorOrder');
                for ii = 1:length(module.Results.History(oo).q)
                    % Cycle through color but periodically reset
                    % when the maximum number of colors is reached.
                    jj = max(mod(ii, size(colorOrder,1)), 1);
                    uq_plot(...
                        module.Results.History(oo).X{ii}(:,1),...
                        module.Results.History(oo).X{ii}(:,2),...
                        'o',...
                        'MarkerSize', 3,...
                        'MarkerFaceColor', colorOrder(jj,:),...
                        'Color', colorOrder(jj,:))
                end
                hold off
                % Set axes labels and title
                xlabel('$\mathrm x_1$') 
                ylabel('$\mathrm x_2$')
                title('\textrm SubsetSim - Samples in each subset')
                
            % Plot the subsets in a histogram for 1-dimensional problems
            case 1
                uq_figure
                X = cell2mat(module.Results.History(oo).X);
                uq_histogram(X)
                xlabel('$\mathrm x$')
                title('SubsetSim - Samples in each subset')
                
            otherwise
                %plot nothing for more than 2 dimensions
                
        end % switch dimensions
        
    end % for oo
    
end %if model evaluations

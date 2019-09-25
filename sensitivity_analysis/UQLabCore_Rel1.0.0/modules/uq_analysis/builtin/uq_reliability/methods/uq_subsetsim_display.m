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
                Colors = [1 0 1
                    0 0 1
                    1 0 0
                    0 1 0
                    1 1 0
                    0 0 0];
                uq_figure; hold on
                for ii = 1:length(module.Results.History(oo).q)
                    uq_plot(module.Results.History(oo).X{ii}(:,1), module.Results.History(oo).X{ii}(:,2), 'o', 'Color', Colors(mod(ii,6)+1,:), 'MarkerSize', 3)
                end; hold off, box on; grid on
                uq_setInterpreters(gca)
                xlabel('$\mathrm x_1$') 
                ylabel('$\mathrm x_2$')
                title('\textrm SubsetSim - Samples in each subset')
                
            % Plot the subsets in a histogram for 1-dimensional problems
            case 1
                uq_figure;
                X = cell2mat( module.Results.History(oo).X );
                uq_histogram(X);
                uq_setInterpreters(gca)
                xlabel('$\mathrm x$')
                title('SubsetSim - Samples in each subset')
                box on; grid on
                
            otherwise
                %plot nothing for more than 2 dimensions
                
        end % switch dimensions
        
    end % for oo
    
end %if model evaluations

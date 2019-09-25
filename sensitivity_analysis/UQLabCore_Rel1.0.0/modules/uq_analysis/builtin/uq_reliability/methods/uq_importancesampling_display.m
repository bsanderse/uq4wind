function varargout = uq_importancesampling_display( module, idx, varargin )
% UQ_IMPORTANCESAMPLING_DISPLAY visualizes importance sampling analysis and 
% its results
%
% See also: UQ_SR_IMPORTANCE_SAMPLING, UQ_DISPLAY_UQ_RELIABILITY

Results = module.Results;

% for each output index
for oo = idx
    
    %% plot the convergence curve for the failure probability estimate
    iter = length(Results.History(oo).Pf);
    if iter ~=1
        N = (1:iter)*module.Internal.Simulation.BatchSize;
        uq_figure
        hold on
        grid on
        f1 = fill([N, fliplr(N)], [Results.History(oo).Pf+Results.History(oo).Conf, fliplr(Results.History(oo).Pf - Results.History(oo).Conf) ],'g');
        set(f1, 'FaceColor', [0.9 0.9 0.9])
        h1 = uq_plot(N, Results.History(oo).Pf, 'LineWidth', 2, 'Color', 'b');
        set(gca, 'xtick', unique(round(get(gca, 'xtick'))))
        l1 = uq_legend([h1,f1], '$\mathrm{P_f}$', 'CI');
        set(l1, 'Interpreter', 'latex')
        uq_setInterpreters(gca)
        xlabel('$\mathrm N$')
        ylabel('$\mathrm{P_f}$')
        title('IS - Convergence', 'interpreter', 'latex')
        box on
        xlim([N(1), N(end)]);
    end
    
    %% display design point, cloud of points in 2 dimensions
    switch length(module.Internal.Input.Marginals)
        case 2
            
            uq_figure()
            hold on
            grid on
            
            %plot the cloud of failed and save importance sampling samples
            if module.Internal.SaveEvaluations
                USamples = Results.History(oo).U;
                LSF = Results.History(oo).G;
                a1 = uq_plot(USamples(LSF<=0,1), USamples(LSF<=0,2), 'or', 'markersize', 3);
                a2 = uq_plot(USamples(LSF>0, 1), USamples(LSF>0, 2), 'og', 'markersize', 3);
            end
            
            % Plot the algorithm steps
            UstarValues = Results.FORM.History(oo).U; % Form steps
            h1 = uq_plot(UstarValues(:,1),UstarValues(:,2),'r->','MarkerFaceColor','r', 'LineWidth', 2);
            
            % Highlight in black the starting point and in green the ending point:
            uq_plot(UstarValues(1,1),UstarValues(1,2),'>k');
            uq_plot(UstarValues(end,1),UstarValues(end,2),'>g');
            
            %plot the FORM limit state surface
            axis equal
            h2 = uq_plot(UstarValues(end,1)+[UstarValues(end,2) -UstarValues(end,2)] ,...
                UstarValues(end,2)+[-UstarValues(end,1), +UstarValues(end,1) ] , 'k', 'LineWidth', 2);
            
            %axis and stuff
            uq_setInterpreters(gca)
            title('IS - FORM design point and failure plane');
            xlabel('$\mathrm u_1$');
            ylabel('$\mathrm u_2$');
            
            box on
            if module.Internal.SaveEvaluations
                l1 = uq_legend([h1, h2, a1, a2], 'FORM iterations', 'FORM limit state surface', '$\mathrm g(X)\leq 0$', '$\mathrm g(X)>0$');
            else
                l1 = uq_legend([h1, h2], 'FORM iterations', 'FORM limit state surface');
            end
            set(l1, 'interpreter', 'latex')
            set(l1, 'Location', 'best')
            
    end
    
end
function varargout = uq_mc_display( module, idx, varargin )
% UQ_MC_DISPLAY visualizes analysis and results of the  Monte Carlo
% simulation
% 
% See also: UQ_DISPLAY_UQ_RELIABILITY

Results = module.Results;

%for each response quantity
for oo = idx
    %only plot when more than one iteration has been done in MCS
    iter = length(Results.History(oo).Pf);
    if iter ~=1
        %% plot the convergence curve for the failure probability estimate 
        N = (1:iter)*module.Internal.Simulation.BatchSize;
        uq_figure
        hold on
        grid on
        f1 = fill([N, fliplr(N)],...
            [Results.History(oo).Pf+Results.History(oo).Conf;...
            flipud(Results.History(oo).Pf - Results.History(oo).Conf) ],...
            'g');
        set(f1, 'FaceColor', [0.9 0.9 0.9])
        h1 = uq_plot(N, Results.History(oo).Pf, 'LineWidth', 2, 'Color', 'b');
        l1 = uq_legend([h1,f1], '$\mathrm{P_f}$', 'CI');
        set(l1, 'Interpreter', 'latex')
        uq_setInterpreters(gca)
        xlabel('$\mathrm{N}$')
        ylabel('$\mathrm{P_f}$')
        title('MCS - Convergence of $\mathrm{P_f}$')
        box on
        xlim([N(1), N(end)])
        
        %% plot the convergence curve for the reliability index 
        uq_figure
        hold on
        grid on
        f2 = fill([N, fliplr(N)],  [-icdf('normal',Results.History(oo).Pf+Results.History(oo).Conf, 0, 1);...
            -icdf('normal', flipud(Results.History(oo).Pf - Results.History(oo).Conf), 0, 1) ],'g');
        set(f2, 'FaceColor', [0.9 0.9 0.9])
        h2 = uq_plot(N, -icdf('normal', Results.History(oo).Pf, 0, 1), 'LineWidth', 2, 'Color', 'b');
        l2 = uq_legend([h2,f2], '$\mathrm{\beta_{MC}}$', 'CI');
        set(l2, 'Interpreter', 'latex')
        uq_setInterpreters(gca)
        xlabel('$\mathrm{N}$')
        ylabel('$\mathrm{\beta_{MC}}$')
        title('MCS - Convergence of $\mathrm{\beta_{MC}}$')        
        box on
        xlim([N(1), N(end)])
    end
    
    switch length(module.Internal.Input.Marginals)
        case 2
            if module.Internal.SaveEvaluations
                
                uq_figure()
                grid on
                
                nplot = min(size(module.Results.History(end).G,1),1e4);
                
                LSF = module.Results.History(end).G(1:nplot,oo);
                XSamples = module.Results.History(end).X(1:nplot,:);
                
                a1 = uq_plot(XSamples(LSF<=0,1), XSamples(LSF<=0,2),...
                    'or', 'MarkerSize', 3);
                hold on
                a2 = uq_plot(XSamples(LSF>0, 1), XSamples(LSF>0, 2),...
                    'og', 'MarkerSize', 3);
                hold off
                
                uq_setInterpreters(gca)
                title('MCS - Samples');
                xlabel('$\mathrm{x_1}$')
                ylabel('$\mathrm{x_2}$')
                
                box on
                if isempty(a1)
                    l1 = uq_legend(a1,'$\mathrm{g(X)\leq 0}$');
                elseif isempty(a2)
                    l1 = uq_legend(a2,'$\mathrm{g(X)>0}$');
                else
                    l1 = uq_legend([a1,a2],...
                        '$\mathrm{g(X)\leq 0}$', '$\mathrm{g(X)>0}$');
                end
                set(l1, 'Interpreter', 'latex', 'Location', 'best')
                
            end
    end
            
end
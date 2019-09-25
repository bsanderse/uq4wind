function varargout = uq_sorm_display( module, idx, varargin ) 
% UQ_SORM_DISPLAY visualizes the analysis and results of SORM and FORM
% analyses
%
% See also: UQ_DISPLAY_UQ_RELIABILITY

Results = module.Results(end);

% for each idx
for oo = idx

%% display the evolution of the reliability index
uq_figure()
hold on
grid on
uq_plot(Results.History(oo).BetaHL, 'bs-', 'LineWidth', 2)
xlim([1 Results.Iterations(oo)])
uq_setInterpreters(gca)
xlabel('number of iterations')
ylabel('$\mathrm{\beta_{HL}}$')
title('FORM - Convergence')
set(gca, 'xtick', unique(round(get(gca, 'xtick'))))
box on

%% display iterations, design point and FORM plane
switch length(module.Internal.Input.Marginals)
    case 2
        % Plot the algorithm steps
        uq_figure()
        hold on
        grid on
        UstarValues = Results.History(oo).U; % Form steps
        h1 = uq_plot(UstarValues(:,1),UstarValues(:,2),'r->','MarkerFaceColor','r', 'LineWidth', 2);
        
        % Highlight in black the starting point and in green the ending point:
        uq_plot(UstarValues(1,1),UstarValues(1,2),'>k');
        uq_plot(UstarValues(end,1),UstarValues(end,2),'>g');
        
        %plot the FORM limit state surface
        h2 = uq_plot(UstarValues(end,1)+[UstarValues(end,2) -UstarValues(end,2)] ,...
            UstarValues(end,2)+[-UstarValues(end,1), +UstarValues(end,1) ], 'k', 'LineWidth', 2);
        
        %axis etc
        xlimits = get(gca, 'XLim');
        axis equal
        xlim(xlimits)
        uq_setInterpreters(gca)
        title('FORM - Design point, failure plane');
        xlabel('$\mathrm u_1$');
        ylabel('$\mathrm u_2$');
        box on
        l1 = uq_legend([h1, h2],...
            'Iterations', 'FORM limit state surface', 'Location', 'best');
        set(l1, 'Interpreter', 'latex')
    otherwise
        % visualization so far only for 2-dimensional input vectors
end

end 
function varargout = uq_akmcs_display(module,idx,varargin)
% UQ_AKMCS_DISPLAY visualizes AK-MCS analysis and its results
%
% See also: UQ_AKMCS, UQ_DISPLAY_UQ_RELIABILITY

% Retrieve the results
Results = module.Results;
MetaModel = module.Internal.AKMCS.MetaModel;
myKrig = Results.(MetaModel);
%% for each output IDX
for oo = idx
    
    %plot the convergence curve of the failure probability estimate
    uq_figure; hold on
    Nstart = Results.History(oo).NInit;
    d = fill([Results.History(oo).NSamples'+Nstart; flipud(Results.History(oo).NSamples'+Nstart)],...
        [Results.History(oo).PfLower'; flipud(Results.History(oo).PfUpper')], 'g');
    set(d, 'FaceColor', [0.9 0.9 0.9])
    a = uq_plot(Results.History(oo).NSamples+Nstart, Results.History(oo).Pf, 'b', 'LineWidth', 2);
    ylimits = get(gca,'ylim');
    xlimits = get(gca,'xlim');
    set(gca,'xlim', [Nstart xlimits(2)])
    c = errorbar((Results.History(oo).NSamples(end)+Nstart), Results.History(oo).Pf(end), ...
        Results.CoV(oo) * norminv(1-module.Internal.Simulation.Alpha/2,0,1) * Results.History(oo).Pf(end), ...
        'k', 'LineWidth', 2);
    l = uq_legend([a,d,c], '$\mathrm{P_f}$', '$\mathrm{P_f^+, P_f^-}$', 'CI MCS');
    set(l, 'interpreter', 'latex')
    uq_setInterpreters(gca)
    % end
    xlabel('Number of samples')
    ylabel('$\mathrm{P_f}$')
    title('AK-MCS - Convergence')

    %% for the 2-dimensional case, plot the safe and failed samples of the
    %experimental design
    if module.Internal.SaveEvaluations
        switch length(module.Internal.Input.nonConst)
            case 2
                uq_figure
                
                hold on
                X = Results.(MetaModel)(oo).ExpDesign.X;
                G = Results.(MetaModel)(oo).ExpDesign.Y;
                a = uq_plot(X(G<=0,1), X(G<=0,2), 'rs', 'LineWidth', 2);
                b = uq_plot(X(G>0, 1), X(G>0, 2), 'g+', 'LineWidth', 2);
                
                
                minX = min(Results.History(oo).MCSample);
                maxX = max(Results.History(oo).MCSample);
                [xx, yy] = meshgrid(linspace(minX(1), maxX(1), 200), linspace(minX(2), maxX(2), 200));
                zz = reshape(uq_evalModel(myKrig(oo), [xx(:), yy(:)]),size(xx));
                [~, cp] = contour(xx,yy,zz, [0 0], 'r', 'linewidth', 1);
                
                % Format the figure
                labels = {'$\mathrm{g(X)\leq 0}$', '$\mathrm{g(X)>0}$', '$\mathrm{g(X)=0}$'};
                pp = [a,b,cp];
                l = uq_legend(pp,labels([~isempty(a) ~isempty(b) ~isempty(c)]) );
                set(l, 'interpreter', 'latex')
                uq_setInterpreters(gca)
                box on
                xlabel('$\mathrm{x_1}$')
                ylabel('$\mathrm{x_2}$')
                title('AK-MCS - Experimental design')
        end 
    end 
    
end
% ================================== PLOTTING and POSTPROCESSING =========================

if (~exist('nout','var'))
    nout = 1;
end


%% Convergence of Sobol indices

% 
% % loop over number of quantities of interest (length of output vector)
% for q=1:nout
%     
%     figure
%     cmap = get(gca,'ColorOrder');
%     
%     if (find(strcmp(methods,'MC')))
%         if(size(AVG_Sobol_MC_Total,2)>1)
%             semilogx(NsamplesMC', AVG_Sobol_MC_Total(:, 1:end-1, q), 'x-','Linewidth', 2, 'Color', cmap(1,:), 'HandleVisibility','off');
%             hold on
%         end
%         semilogx(NsamplesMC', AVG_Sobol_MC_Total(:, end, q), 'x-','Linewidth', 2, 'Color', cmap(1,:));
%     end
%     if (find(strcmp(methods,'PCE_Quad')))
%         if(size(Sobol_Quad_Total,2)>1)
%             semilogx(NsamplesQuad', Sobol_Quad_Total(:, 1:end-1, q), 's-','Linewidth', 2,'Color', cmap(2,:), 'HandleVisibility','off');
%             hold on
%         end
%         semilogx(NsamplesQuad', Sobol_Quad_Total(:, end, q), 's-','Linewidth', 2,'Color', cmap(2,:));
%     end
%     if (find(strcmp(methods,'PCE_OLS')))
%         if(size(AVG_Sobol_OLS_Total,2)>1)
%             semilogx(NsamplesOLS, AVG_Sobol_OLS_Total(:, 1:end-1, q), 'o-','Linewidth', 2,'Color', cmap(3,:), 'HandleVisibility','off');
%             hold on
%         end
%         semilogx(NsamplesOLS, AVG_Sobol_OLS_Total(:, end, q), 'o-','Linewidth', 2,'Color', cmap(3,:));
%     end
%     if (find(strcmp(methods,'PCE_LARS')))
%         if(size(AVG_Sobol_LARS_Total,2)>1)
%             semilogx(NsamplesLARS, AVG_Sobol_LARS_Total(:, 1:end-1, q), 'd-','Linewidth', 2,'Color', cmap(4,:), 'HandleVisibility','off');
%             hold on
%         end
%         semilogx(NsamplesLARS, AVG_Sobol_LARS_Total(:,end, q), 'd-','Linewidth', 2,'Color', cmap(4,:));
%     end
%     xlabel('N') % Add proper labelling and a legend
%     legend(methods, 'Interpreter', 'none')
%     ylabel('Total index');
%     grid on;
%     title(strcat('Convergence of Sobol indices for output ',num2str(q)))
%     
% end


%% bar chart of Sobol indices

% corresponding to largest number of samples

for q=1:nout
    
    figure
    cmap = get(gca,'ColorOrder');

    hold on
    
    n_methods = length(methods);
    bar_width = 0.5/n_methods;
    bar_vec   = 1:n_methods;
    coords    = (bar_vec - mean(bar_vec))*bar_width;
    k         = 1;
    
    if (find(strcmp(methods,'MC')))
%         if(length(NsamplesMC)==1)
%             uq_bar((1:ndim)+ coords(k), AVG_Sobol_MC_Total(:,end), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
%         else
            uq_bar((1:ndim)+ coords(k), AVG_Sobol_MC_Total(end,:,q), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
%         end
        k = k+1;
    end
    
    if (find(strcmp(methods,'PCE_Quad')))
        uq_bar((1:ndim)+ coords(k), Sobol_Quad_Total(end,:,q), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
        k = k+1;
    end
    
    if (find(strcmp(methods,'PCE_OLS')))
%         if(length(NsamplesOLS)==1)
%             uq_bar((1:ndim)+ coords(k), AVG_Sobol_OLS_Total(:,end), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
%         else
            uq_bar((1:ndim)+ coords(k), AVG_Sobol_OLS_Total(end,:,q), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
%         end
        k = k+1;
    end
    
    if (find(strcmp(methods,'PCE_LARS')))
%         if(length(NsamplesLARS)==1)
%             uq_bar((1:ndim)+ coords(k), AVG_Sobol_LARS_Total(:,end,q), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
%         else
            uq_bar((1:ndim)+ coords(k), AVG_Sobol_LARS_Total(end,:,q), bar_width, 'FaceColor', cmap(k,:), 'EdgeColor', 'none')
%         end
        k = k+1;
    end
    
    legend(methods, 'Interpreter', 'none')
    ylabel('Total order Sobol index');
    ylim([0 1])
    xticks(1:ndim)
    for i =1:ndim
        label_names{i} = Input.Marginals(i).Name;
    end
    xticklabels(label_names)
    title(strcat('Sobol indices for output ',num2str(q)))

end

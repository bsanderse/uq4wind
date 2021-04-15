%% ================================== PLOTTING and POSTPROCESSING =========================
% uq_figure('Position', [50 50 500 400]);
% myColors = uq_cmap(2);

%% mean
if (compare_mean == 1)
    
    figure
    if (find(strcmp(methods,'MC')))
        loglog(NsamplesMC, err_mean_MC, 'x-','Linewidth', 2); %, 'Color', myColors(1,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_Quad')))
        loglog(NsamplesQuad, err_mean_Quad, 's-','Linewidth', 2);%,'Color', myColors(2,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_OLS')))
        loglog(NsamplesOLS, err_mean_OLS, 'o-','Linewidth', 2); %,'Color', 'y');
        hold on
    end
    if (find(strcmp(methods,'PCE_LARS')))
        loglog(NsamplesLARS, err_mean_LARS, 'd-','Linewidth', 2); %,'Color', 'r');
        hold on
    end
    xlabel('N') % Add proper labelling and a legend
    legend({'Monte Carlo','PCE - quadrature','PCE - least squares','PCE - LARS'})
    grid on;
    title('Error in mean')
end


%% standard deviation
if (compare_std == 1)
    
    figure
    if (find(strcmp(methods,'MC')))
        loglog(NsamplesMC, err_std_MC, 'x-','Linewidth', 2); %, 'Color', myColors(1,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_Quad')))
        loglog(NsamplesQuad, err_std_Quad, 's-','Linewidth', 2);%,'Color', myColors(2,:));
        hold on
    end
    if (find(strcmp(methods,'PCE_OLS')))
        loglog(NsamplesOLS, err_std_OLS, 'o-','Linewidth', 2); %,'Color', 'y');
        hold on
    end
    if (find(strcmp(methods,'PCE_LARS')))
        loglog(NsamplesLARS, err_std_LARS, 'd-','Linewidth', 2); %,'Color', 'r');
        hold on
    end
    xlabel('N') % Add proper labelling and a legend
    legend({'Monte Carlo','PCE - quadrature','PCE - least squares','PCE - LARS'})
    grid on;
    title('Error in standard deviation')
    
end


%% Convergence of Sobol indices
figure
cmap = get(gca,'ColorOrder');


if (find(strcmp(methods,'MC')))
%     if(size(AVG_Sobol_MC_Total,2)>1)
%         semilogx(Sobol_MC_Nsamples', AVG_Sobol_MC_Total(:, 1:end-1), 'x-','Linewidth', 2, 'Color', cmap(1,:), 'HandleVisibility','off');
%         hold on
%     end
    semilogx(NsamplesMC', AVG_Sobol_MC_Total(:, end), 'x-','Linewidth', 2, 'Color', cmap(1,:));
end
if (find(strcmp(methods,'PCE_Quad')))
%     if(size(AVG_Sobol_Quad_Total,2)>1)
%         semilogx(NsamplesQuad', AVG_Sobol_Quad_Total(:, 1:end-1), 's-','Linewidth', 2,'Color', cmap(2,:), 'HandleVisibility','off');
%         hold on
%     end
    semilogx(NsamplesQuad', AVG_Sobol_Quad_Total(:, end), 's-','Linewidth', 2,'Color', cmap(2,:));
end
if (find(strcmp(methods,'PCE_OLS')))
%     if(size(AVG_Sobol_OLS_Total,2)>1)
%         semilogx(NsamplesOLS, AVG_Sobol_OLS_Total(:, 1:end-1), 'o-','Linewidth', 2,'Color', cmap(3,:), 'HandleVisibility','off');
%         hold on
%     end
    semilogx(NsamplesOLS, AVG_Sobol_OLS_Total(:, end), 'o-','Linewidth', 2,'Color', cmap(3,:));
end
if (find(strcmp(methods,'PCE_LARS')))
%     if(size(AVG_Sobol_LARS_Total,2)>1)
%         semilogx(NsamplesLARS, AVG_Sobol_LARS_Total(:, 1:end-1), 'd-','Linewidth', 2,'Color', cmap(4,:), 'HandleVisibility','off');
%         hold on
%     end
    semilogx(NsamplesLARS, AVG_Sobol_LARS_Total(:,end), 'd-','Linewidth', 2,'Color', cmap(4,:));
end
xlabel('N') % Add proper labelling and a legend
legend(methods, 'Interpreter', 'none')
ylabel('Total index');
grid on;
title('Comparison of Sobol indices')

%% bar chart of Sobol indices
figure
cmap = get(gca,'ColorOrder');

hold on


m_plot = 1;
n_plot = N_out; % number of columns = number of radial sections

QoI_names = {'induced velocity'};
titles = {'section 1', 'section 2', 'section 3', 'section 4', 'section 5'};

for q=1:nout
    
    [i_plot,j_plot] = ind2sub([m_plot n_plot],q);
    % note that subplot index does not correspond to ind2sub indexing
    % therefore, get q_plot by reversing the indexing:
    q_plot = sub2ind([n_plot m_plot],j_plot,i_plot);
    subplot(m_plot,n_plot,q_plot);
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
    
    %     legend(methods, 'Interpreter', 'none')
    %     ylabel('Total order Sobol index');
    ylim([0 1])
    xticks(1:ndim)
    for i =1:ndim
        label_names{i} = Input.Marginals(i).Name;
    end
    xticklabels(label_names)
    if (j_plot==1) % first column
        ylabel(QoI_names{i_plot});
    end
    if (i_plot==1) % first row
        title(titles{j_plot});
    end
    
end


%% plot the QoI at each radial section
% assuming LARS is used
mySurrogateModel = myPCE_LARS;

% Y has size Nsamples * radial sections
QoI_sim = myPCE_LARS.ExpDesign.Y;

% get unperturbed model result from surrogate if not available yet
% if (test_run == 0) % test run has not yet been performed
%     warning('unperturbed AeroModule results are obtained from the surrogate model');
QoI_unperturbed = uq_evalModel(myModel,X_unperturbed);
%     QoI_unperturbed = Y_unperturbed(:,j:n_coeffs:end);
% elseif (test_run == 1)
%     QoI_unperturbed = Y_unperturbed(:,j:n_coeffs:end);
% end

figure
t  = lines; %get(gca,'ColorOrder');
h2 = plot(r_rel,QoI_sim,'o','markersize',12,'color',t(2,:),'Linewidth', 2.5);
hold on
h3 = plot(r_rel,QoI_unperturbed,'s','markersize',16,'color',t(3,:),'Linewidth', 2.5);
grid on
xlabel('r [m]');
%         ylabel('Fn [N/m]');
%     legend([h2 h3],'Perturbations','Unperturbed');
%     title(['Model vs. data for run ' num2str(select_runs(i)) ' and Fourier coefficient ' num2str(k)]);




%% plot the polynomial response surface for PCE quadrature-based

if (find(strcmp(methods,'PCE_Quad')))
    
    n_inputs = length(Input.Marginals);
    
    if (n_inputs == 1)
        Xsamples  = myPCE_Quad.ExpDesign.X;
        Ysamples  = myPCE_Quad.ExpDesign.Y;
        
        N = 100;
        domain = getMarginalBounds(myInput.Marginals(1));
        X      = linspace(domain(1),domain(2),N)';
        
        Y     = uq_evalModel(myModel,X);
        Y_PCE = uq_evalModel(myPCE_Quad,X); % or myPCE_OLS, myPCE_LARS
        
        figure
        plot(X,Y,'-');
        hold on
        plot(X,Y_PCE,'--');
        plot(Xsamples,Ysamples,'s');
        
        
    elseif (n_inputs >= 2)
        
        Xsamples  = myPCE_Quad.ExpDesign.X;
        Ysamples  = myPCE_Quad.ExpDesign.Y;
        
        % plot two parameters
        N1 = 100;
        N2 = 100;
        
        % select which parameters to plot
        p1 = 1;
        p2 = 2;
        p3 = 3; % we will use mean for p3
        
        % create regular grid of points where surrogate model will be evaluated
        % (this should be cheap)
        domain1 = getMarginalBounds(myInput.Marginals(p1));
        X1      = linspace(domain1(1),domain1(2),N1)';
        domain2 = getMarginalBounds(myInput.Marginals(p2));
        X2      = linspace(domain2(1),domain2(2),N2)';
        
        
        if (n_inputs==2)
            X1new=kron(X1,ones(N2,1));
            X2new=kron(ones(N1,1),X2);
            X = [X1new X2new];
            plotting = 1;
        elseif (n_inputs==3)
            % for n>2, it is difficult to plot something in 3D space
            % for quadrature methods based on tensor formulation,
            % number of samples in one direction = degree+1
            % take 'middle' (mean?) point of samples in X3
            samples1D = 1+DegreesQuad(end);
            start1D  = ceil(samples1D/2);
            X3    = Xsamples(start1D,p3);
            
            X1new = kron(X1,ones(N2,1));
            X2new = kron(ones(N1,1),X2);
            X3new = ones(size(X1new))*X3;
            X = [X1new X2new X3new];
            plotting = 1;
        else
            warning('plotting not correctly implemented for high dimensional spaces');
            plotting = 0;
        end
        
        if (plotting == 1)
            % evaluate original model (can be expensive)
            Y     = uq_evalModel(myModel,X);
            % evaluate surrogate model (should be cheap)
            Y_PCE = uq_evalModel(myPCE_Quad,X); % or myPCE_OLS, myPCE_LARS
            
            figure
            surf(X1,X2,reshape(Y_PCE,N2,N1));
            hold on
            %     surf(X1,X2,reshape(Y,N2,N1));
            if (n_inputs==2)
                plot3(Xsamples(:,1),Xsamples(:,2),Ysamples,'s','MarkerSize',16,'MarkerFaceColor','black');
            elseif (n_inputs==3)
                plot3(Xsamples(start1D:samples1D:end,1),Xsamples(start1D:samples1D:end,2),Ysamples(start1D:samples1D:end),'s','MarkerSize',16,'MarkerFaceColor','black');
            end
            title('polynomial response surface');
        end
        
    end
    
end
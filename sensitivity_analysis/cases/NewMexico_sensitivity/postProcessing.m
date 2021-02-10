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

figure
cmap = get(gca,'ColorOrder');

hold on

m_plot = 2; % number of coefficients used for QoI
n_plot = 5; % number of columns = number of radial sections

titles = {'section 1', 'section 2', 'section 3', 'section 4', 'section 5'};
% QoI_names   = {'mean','amplitude 1','angle 1'};
QoI_names   = {'cosine','sine'};

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
%     title(strcat('Sobol indices for output ',num2str(q)))

end

%% convergence of LOO and modLOO

figure
for k=1:nout

    if (length(NsamplesLARS)==1)
        semilogy(NsamplesLARS,AVG_LOO_LARS(k),'s-','LineWidth',2);
    else
        semilogy(NsamplesLARS,AVG_LOO_LARS(:,k),'s-','LineWidth',2);
    end
    hold on
end
grid on
switch fourier_type
    
    case 'real_imag'
        legend('Cosine Section 1','Cosine Section 2','Cosine Section 3','Cosine Section 4','Cosine Section 5',...
               'Sine Section 1','Sine Section 2','Sine Section 3','Sine Section 4','Sine Section 5');
    case 'amp_phase'
         legend('AM Section 1','AM Section 2','AM Section 3','AM Section 4','AM Section 5',...
                'PH Section 1','PH Section 2','PH Section 3','PH Section 4','PH Section 5');
    otherwise
        warning('wrong fourier_type selected');
end
xlabel('Number of samples (=number of model runs)'); 
ylabel('LOO error');


figure
for k=1:nout

    if (length(NsamplesLARS)==1)
        semilogy(NsamplesLARS,AVG_modLOO_LARS(k),'s-','LineWidth',2);
    else
        semilogy(NsamplesLARS,AVG_modLOO_LARS(:,k),'s-','LineWidth',2);
    end
    hold on
end
grid on
switch fourier_type
    
    case 'real_imag'
        legend('Cosine Section 1','Cosine Section 2','Cosine Section 3','Cosine Section 4','Cosine Section 5',...
               'Sine Section 1','Sine Section 2','Sine Section 3','Sine Section 4','Sine Section 5');
    case 'amp_phase'
         legend('AM Section 1','AM Section 2','AM Section 3','AM Section 4','AM Section 5',...
                'PH Section 1','PH Section 2','PH Section 3','PH Section 4','PH Section 5');
    otherwise
        warning('wrong fourier_type selected');
end
xlabel('Number of samples (=number of model runs)'); 
ylabel('Modified LOO error');

%% plot the actual response surface

% this is easy if we have one uncertain parameter:
if (nunc == 1)
    
    
    if (find(strcmp(methods,'PCE_LARS')))        
            
            % points where the code has been evaluated:
            % note that Xsamples includes the constants
            Xsamples  = myPCE_LARS.ExpDesign.X;
            Ysamples  = myPCE_LARS.ExpDesign.Y;            
            
%             Xsamples  = myPCE_LARS.ExpDesign.X;
%             Ysamples  = unwrap(myPCE_LARS.ExpDesign.Y);
%             
%             % adapt the PCE by unwrapping the angle dependent part
%             metamodel_LARS_unwrapped = metamodelLARS;
%             %change the experimental design
%             metamodel_LARS_unwrapped.ExpDesign.X = myPCE_LARS.ExpDesign.X;
%             metamodel_LARS_unwrapped.ExpDesign.Y = myPCE_LARS.ExpDesign.Y;
%             metamodel_LARS_unwrapped.ExpDesign.Y(:,2:2:end) = Ysamples(:,2:2:end);
%             metamodel_LARS_unwrapped.ExpDesign.Sampling='user';
%             myPCE_LARS_unwrapped  = uq_createModel(metamodel_LARS_unwrapped);

            
            % evaluate surrogate model at many points:
            Ntest = 100;
            domain = getMarginalBounds(myInput.Marginals(1));
            
            X_PCE = linspace(domain(1),domain(2),Ntest)';   
            % add constants: 
            X_PCE_full = [X_PCE repmat(Xsamples(1,nunc+1:end),Ntest,1)];

            % evaluate the PCE model at many points
            Y_PCE = uq_evalModel(myPCE_LARS,X_PCE_full); 
            
            % assume that amplitude and phase of 1 Fourier mode are
            % considered
            % amplitude response at different sections
            figure
            cmap = get(gca,'ColorOrder');

            for k=1:n_r_index
                plot(X_PCE(:,1:nunc),Y_PCE(:,2*k-1),'-','Color',cmap(k,:));
                hold on
                plot(Xsamples(:,1:nunc),Ysamples(:,2*k-1),'s','Color',cmap(k,:));
            end
            xlabel(myInput.Marginals(1).Name);
%             ylabel('amplitude first Fourier mode')
            ylabel('cosine coefficient first Fourier mode')
            legend('section 1 - surrogate','section 1 - AeroModule run',...
                'section 2 - surrogate','section 2 - AeroModule run',...
                'section 3 - surrogate','section 3 - AeroModule run',...
                'section 4 - surrogate','section 4 - AeroModule run',...
                'section 5 - surrogate','section 5 - AeroModule run');
            grid on
            
            % phase response at different sections
            figure
            cmap = get(gca,'ColorOrder');
            for k=1:n_r_index
                plot(X_PCE(:,1:nunc),Y_PCE(:,2*k),'-','Color',cmap(k,:));
                hold on
                plot(Xsamples(:,1:nunc),Ysamples(:,2*k),'s','Color',cmap(k,:));            
            end
            xlabel(myInput.Marginals(1).Name);
%             ylabel('phase shift first Fourier mode') 
            ylabel('sine coefficient first Fourier mode')

            legend('section 1 - surrogate','section 1 - AeroModule run',...
                'section 2 - surrogate','section 2 - AeroModule run',...
                'section 3 - surrogate','section 3 - AeroModule run',...
                'section 4 - surrogate','section 4 - AeroModule run',...
                'section 5 - surrogate','section 5 - AeroModule run');
            grid on
    else            
            disp('response surface plotting not implemented for this method');
    end
    
    
else
    
    disp('response surface plotting not implemented for more than 1 uncertain parameter');
    
end

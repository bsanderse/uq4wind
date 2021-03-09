clearvars
close all

folder_list = {'10regularsamples_mean','50regularsamples_mean','100regularsamples_mean','200regularsamples_mean'};

% choose which variable to plot (give index):
X_choice  = 5;
plot_prior = true;

for i=1:length(folder_list)
    
    disp(folder_list{i});
    data = load(fullfile(folder_list{i},'workspace.mat'));
    BayesianAnalysis = data.BayesianAnalysis;
    Prior = data.BayesianAnalysis.PriorDist;
    
    %% process the Bayesian Analysis
    uq_postProcessInversion(BayesianAnalysis,'burnIn',0.5,'posteriorPredictive',0)
    
    
    %% plot the posterior for a certain variable
    
    % all posterior samples (except the ones disregarded with burn-in)
    % posterior is an array of size Nsteps*Nvars*Nchains
    posterior = BayesianAnalysis.Results.PostProc.PostSample;
    % select the wanted variable, and reshape the rest into a 1D array
    post_X_choice = reshape(squeeze(posterior(:,X_choice,:)),[],1);
    % make density estimate:
    [post_density, x_post_density] = ksdensity(post_X_choice);
    
    if (i==1 && plot_prior)
        % add the prior density:
        switch Prior.Marginals(X_choice).Type
            
            case 'Uniform'
                a = Prior.Marginals(X_choice).Parameters(1);
                b = Prior.Marginals(X_choice).Parameters(2);
                x_prior = linspace(a,b,100);
                prior   = (x_prior>=a & x_prior<=b)/(b-a);
                
                
            case 'Gaussian'
                
                if (Prior.Marginals(X_choice).Bounds(1)) % bounds -> truncated gaussian
                    
                    % truncated gaussian:
                    x_prior = linspace(Prior.Marginals(X_choice).Bounds(1),Prior.Marginals(X_choice).Bounds(2),100);
                    [prior,x_prior] = truncnormpdf(x_prior,...
                        Prior.Marginals(X_choice).Moments(1),Prior.Marginals(X_choice).Moments(2),...
                        Prior.Marginals(X_choice).Bounds(1),Prior.Marginals(X_choice).Bounds(2));
                else
                    mu = Prior.Marginals(X_choice).Moments(1);
                    sigma = Prior.Marginals(X_choice).Moments(2);
                    x_prior = linspace(mu-2*sigma,mu+2*sigma,100);
                    prior   = normpdf(x_prior,mu,sigma);
                end
                
                
            otherwise
                error('prior type not known');
        end
        
    end
    figure(1)
    if (i==1 && plot_prior)
        plot(x_prior,prior,'LineWidth',2);
    end
    hold on
    plot(x_post_density,post_density,'LineWidth',2);
    
    clear data;
    
end

%%
legend('prior','10 measurements','50 measurements','100 measurements','200 measurements');
grid on
xlabel('\Delta C_{l,1}');
ylabel('probability density');
set(gca,'FontSize',12);
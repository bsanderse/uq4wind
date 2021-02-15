%% This file is used to initialize all settings in order to perform the calibration of the AeroModule 
% for the case that experimental data corresponding to different operating
% conditions is available
% For 3D correction calibration please set the desired correction type in
% the respective field below


%% Turbine name and folders

% Name of Matlab file that describes the uncertain parameters used for calibration
% check this file for definition of uncertainties
turbineName = 'NewMexico_calibrate'; 

% folder used to get reference input files
ref_folder      = 'AEROmodule/NM_cal_TWS/reference/';
% folder used to write new input files
current_folder  = 'AEROmodule/NM_cal_TWS/current/';
% name of input file is assumed to be input.txt


%% Forward model description 
% the following settings are valid independent of the operating conditions

% Name of Matlab file representing the model, this is the AeroModule 
% which is called in aero_module.m
Model.mHandle = @aero_module;
% Quantity of interest (see NewMexico_calibrate_readoutput)
QoI = 'Sectional_normal_force';
% choose whether only the time-average of the data is used, or the full dataset
% using 'mean' or 'full'
% for the full dataset, a Fourier analysis is done, and one should choose
% the number of fourier modes below
% note that the 'mean' case can also be run by setting 'full' and then
% choosing n_fourier=1.
QoI_type = 'full';

%% Select Correction type, type = 1 -->Snel, type = 2 --> Chaviaropoulos - Hansen, 0 no external correction,
% for no correction need to set  3d corr flag in input.txt to 1 so Snel is applied

cor_type =0;
%%
% the following settings are only used in case of 'full':
% (currently for 'mean', only 1 revolution is used, and all radial indices)
% number of revolutions to consider (counting from end of time series)
n_rev = 4;
% number of Fourier coefficients to keep (including mean)
% note: we get (n_fourier-1)*2 + 1 coefficients since there is both a real and
% imaginary component (stored as amplitude and phase angle) for each
% frequency
n_fourier = 1;
% radial indices (blade sections) to consider:
r_index = 1:5;
    

%% High-level surrogate model options

% perform test run with Forward Model without uncertainties; 
% this is used as a check and to plot the uncalibrated model in the results
test_run = 1; 

% Switch for Bayesian analysis with the AeroModule or with the surrogate model
Bayes_full = 0; % 0: use and/or set-up surrogate model (PCE); 1: run full model for Bayes (Computationally expensive!)

% If Bayes_full = 0, we need to specify options for loading a surrogate model
Surrogate_model_type = 1; % 0: Uses a stored PCE surrogate model, 1: create surrogate model

% Options for loading a surrogate model
% note that this file should contain all the surrogate models
Surrogate_model_filename = 'StoredSurrogates/NM_cal_TWS/PCE_testrun.mat'; 
   

%% Detailed surrogate model options

% These settings are used if Bayes_full = 0 and Surrogate_model_type = 1
MetaOpts.Type     = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method   = 'LARS'; % Quadrature, OLS, LARS

MetaOpts.ExpDesign.Sampling = 'LHS';
MetaOpts.ExpDesign.NSamples = 100; % number of samples for each surrogate model (each run)
MetaOpts.Degree = 1:4;
MetaOpts.TruncOptions.qNorm = 0.75;   


%% Description of experimental dataset and choice of operating conditions

% NewMexicoData as obtained from Koen Boorsma (TNO)
% folder that contains the experimental data
folder_exp    = fullfile(root_folder,'..','Experimental','NewMexicoData');
% filename with description of datasets
filename_runs = fullfile(folder_exp,'DPN_overview.csv');
% choose the conditions (columns) from the table that are used in the input files
% these should be strings that match the column headers
changing_conditions = {'AIRDENSITY','PITCHANGLE','YAWANGLE','WINDSPEED'}; 
% choose the runs that are to be included in the calibration
% for all runs, set select_runs = 928:957;
select_runs = [432 439 447] ; 

% the position of the sections of the experimental data which are used for
% interpolation of the aeromodule results: see NewMexico_calibrate_readoutput.m
r_exp_data = [0.25 0.35 0.6 0.82 0.92]*2.25;


%% MCMC options
% MCMC parameters
Solver.Type = 'MCMC';
% MCMC algorithms available in UQLab
MCMC_type = 'AIES'; % 'MH': Metropolis-Hastings, 'AM': Adaptive Metropolis
% 'AIES': Affine invariant ensemble, 'HMC': Hamilton Monte Carlo

switch MCMC_type
    
    case 'MH'
        
        Solver.MCMC.Sampler = 'MH';
        Solver.MCMC.Steps = 1e2;
        Solver.MCMC.NChains = 1e2;
        Solver.MCMC.T0 = 1e1;
        
    case 'AM'
        Solver.MCMC.Sampler = 'AM';
        Solver.MCMC.Steps = 1e2;
        Solver.MCMC.NChains = 1e2;
        Solver.MCMC.T0 = 1e1;
        Solver.MCMC.Epsilon = 1e-2;
        
    case 'AIES'
        Solver.MCMC.Sampler = 'AIES';
        Solver.MCMC.Steps = 1e2;
        Solver.MCMC.NChains = 1e1;
        Solver.MCMC.a = 5;
        
    case 'HMC'
        Solver.MCMC.Sampler = 'HMC';
        Solver.MCMC.LeapfrogSteps = 1e3;
        Solver.MCMC.LeapfrogSize = 0.01;
        Solver.MCMC.Mass = 100;
        
    otherwise
        error('wrong MCMC type provided');
        
end


%% start constructing the model structure as used in UQLab required for the calibration

% Get and store uncertain and fixed parameters
% Pass parameters to model via the cell array FixedInputs
% these are the uncertain inputs without the operating conditions
[FixedParameters,UncertainInputs_NoOC] = getParameterAeroModule(turbineName);
   
FixedParameters.root_folder    = root_folder; % working directory
FixedParameters.ref_folder     = ref_folder;  % folder where the AeroModule reference files are located
FixedParameters.current_folder = current_folder;  % folder where the AeroModule adapted files will be located

FixedParameters.QoI            = QoI;
FixedParameters.QoI_type       = QoI_type;
FixedParameters.correction     = cor_type;

switch QoI_type
    case 'full'
        % store parameters in struct
        FixedParameters.n_rev          = n_rev;
        FixedParameters.n_fourier      = n_fourier;
        FixedParameters.r_index        = r_index;
end

P.FixedParameters = FixedParameters;

%% 
if (Bayes_full == 0) % use a PCE surrogate model for calibration
    if (Surrogate_model_type == 0) % load existing surrogate model
        disp(['loading surrogate models from file: ' Surrogate_model_filename]);
        loaded_surrogate_models = load(Surrogate_model_filename);
        % check whether loaded surrogate model is having same features as
        % the uncertainties given in the prior (those we are trying to calibrate)
        if (length(loaded_surrogate_models.mySurrogateModels)~=length(select_runs))
            error('Number of surrogate model loaded from file is different from required number');
        end
    end
end


%% now loop over all selected runs:
% - set up a forward and surrogate model for each operating condition
% - get the experimental data for this operating conditions
% - set up likelihood and prior

% get the experimental conditions into a table
data_runs     = readtable(filename_runs,'PreserveVariableNames',true);
% determine number of experimental runs
% n_runs = size(data_runs,1);
% number of runs chosen
n_runs = length(select_runs);
% number of columns chosen
n_oper = length(changing_conditions);
% total number of columns in table
n_cols = size(data_runs,2); 

% loop over all selected runs
for i = 1:n_runs
    
    %% get the operating conditions for the current run    
    % find the index in the table of the current run
    i_run = find(data_runs.PNT == select_runs(i));
    
    if (isempty(i_run))
        error(['run with ID ' num2str(select_runs(i)) ' not found']);
    end
    
    % loop over the columns of the table and build an object
    % OperatingConditions for this run that contains the selected operating
    % conditions
    k = 1;
    conditions_covered = zeros(n_oper,1);
    for col=1:n_cols
        
        % select the column that matches the selected operating conditions
        id_var = strcmp(data_runs.Properties.VariableNames{col},changing_conditions);
 
        % store in OperatingConditions
        if (any(id_var))
            ind = find(id_var==1);
            OperatingCondition(k).Name = data_runs.Properties.VariableNames{col};
            OperatingCondition(k).Parameters = data_runs(i_run,col).Variables;
            conditions_covered(k) = 1;
            k=k+1;
            
        end
        
    end
    ind_notcovered = find(conditions_covered==0);
    if (~isempty(ind_notcovered))
        for k=1:length(ind_notcovered)
            warning(['Variable ' changing_conditions{k} ' not found in the table']);
        end
    end
    
    % add Operatingconditions to the structure that contains the
    % uncertainties
    UncertainInputs = addOperatingConditions(UncertainInputs_NoOC,OperatingCondition);
    clear OperatingConditions; 
    
    % total dimension of the uncertain inputs vector, this includes any
    % constants
    ndim  = length(UncertainInputs.Marginals);

    
    %% get the experimental dataset for the current run   

    % filename uses the R, P and D columns of the csv file
    filename_exp = strcat('R',num2str(data_runs.RUN(i_run)),'P',num2str(data_runs.POL(i_run)),'D',num2str(data_runs.PNT(i_run)),'_loads.dat');
    full_filename_exp = fullfile(folder_exp,filename_exp);
    % read in the table
    output_raw   = readNewMexico(full_filename_exp);    

    % azimuthal positions
    azi_exp_data = output_raw.Azi;
    
    % Because the model has different discrepancy options at different radial locations,
    % the measurement data is stored in five different data structures:
    % normal forces at five stations
    Fn_exp_data  = table2array(output_raw(:,2:6));
    
    switch QoI_type
        
        case 'full'
            % concatenate all the time-dependent normal forces into one row vector :
            % [ Fn1(t) Fn2(t) Fn3(t) Fn4(t) Fn5(t)]
            dazi  = 10;
%             dt    = dazi/RPM/6;
            azi_exp_int = (0:dazi:350)';
            Fn_exp_int  = spline(azi_exp_data',Fn_exp_data',azi_exp_int)';
            % get the coefficients of the first n_fourier modes
            % the coefficients are ordered according to the PSD
            Fhat        = getFourierCoefficients(Fn_exp_int,n_fourier);
            ind_select  = (2:2:2*(n_fourier-1))';
            
            Fhat_total = [];
            n_r_index = length(r_index);
            n_coeffs  = 2*n_fourier-1; % mean + ampl. mode 1 + angle mode 1 + ampl. mode 2 + angle mode 2 + etc.
            for k = 1:n_r_index
                
                Fhat_mean   = abs(Fhat(1,r_index(k))); 
                Fhat_new    = Fhat(ind_select,r_index(k));  
                Fhat_total  = horzcat(Fhat_total,[Fhat_mean 2*abs(Fhat_new).' angle(Fhat_new).']);
                
            end
            Data(i).y    = Fhat_total;            
            Data(i).Name = 'Normal force';
                        
            % describe the model ID and output ID corresponding to this
            % data set
            % in this case, Data(i) corresponds to the i-th model (i-th operating condition)
            % to the Model ID map is simple
            % the output ID maps the entry in Data(i) to the entries in the model output Y
            % (as computed in aero_module.m)
            % in our case, Y is constructed such that it is ordered in
            % the same way as Data.y
            n_output      = n_coeffs * n_r_index;
            Data(i).MOMap = [ i*ones(1,n_output); % Model IDs ...
                              1:n_output]; % Output IDs
            
        case 'mean'
            
            % row vector
            n_coeffs = 1;
            Data(i).y = mean(Fn_exp_data);
            Data(i).Name = 'Normal force';
            Data(i).MOMap = [ i i i i i; % Model ID ...
                              1 2 3 4 5]; % Output ID
            
    end  
    % prevent output_raw from being reused
    clear output_raw;     
    
    
    %% set-up the forward model

    % we need some of the experimental locations to interpolate the output
    P.FixedParameters.r_exp          = r_exp_data;
    P.FixedParameters.azi_exp        = azi_exp_data;    
    % set uncertain inputs to the model
    P.UncertainInputs = UncertainInputs;
    Model.Parameters  = P;
    %  note that the model output is (in general) a vector, containing e.g. the
    %  force at different radial sections; this does not require any additional
    %  specification in UQLab
    %  the model is however not vectorized in the sense that we cannot give all
    %  the possible parameter values at once as input to AeroModule, but
    %  instead we need to do this sequentially
    Model.isVectorized = false;    
    
    Model.Name     = ['Normal force ' num2str(select_runs(i))];
    
    % create the forward model
    myForwardModel = uq_createModel(Model);
    
    %% Prior
    % Set the Prior equal to the Input
    Prior   = UncertainInputs;
    myPrior = uq_createInput(Prior);
    

    % do a test run with the forward model at unperturbed settings
    ndim = length(myPrior.Marginals);
    % set unperturbed vector:
    for k=1:ndim
        % we take the mean of each parameter as the unperturbed
        % condition
        X_unperturbed(1,k) = myPrior.Marginals(k).Moments(1);
    end    
    if (exist('test_run','var'))
        if (test_run == 1)
            disp('Performing test run at unperturbed (mean value) settings');
            % store output of current model i
            Y_unperturbed(i,:) = uq_evalModel(myForwardModel,X_unperturbed);
        end
    end        
    
    %% set-up surrogate model (if needed)
    
    if (Bayes_full == 0) % use a PCE surrogate model for calibration
        if (Surrogate_model_type == 0)
            % check whether loaded surrogate model is having same features as
            % the uncertainties given in the prior (those we are trying to calibrate)
            if (~isequaln(myPrior.Marginals,loaded_surrogate_models.mySurrogateModels(i).Model.Options.Input.Marginals))
                error('Marginals specified for Prior do not correspond with marginals used in surrogate model');
            end
            
        elseif (Surrogate_model_type == 1)
            disp(['creating surrogate model for run ' num2str(select_runs(i))]);
            % use prior also as input uncertainties for surrogate model
            MetaOpts.Input     = myPrior;
            MetaOpts.FullModel = myForwardModel;
            MetaOpts.Name      = strcat(Model.Name,'-surrogate');

            % create a surrogate model for each run (so for each operating
            % condition)
            mySurrogateModels(i).Model   = uq_createModel(MetaOpts);    
            % each model uses the entire parameter vector:
            mySurrogateModels(i).PMap    = 1:ndim;
            
        end
    else % do Bayes with full model (not recommended)
        
        % add the current forward model into the myForwardModels structure
        myForwardModels(i).Model = myForwardModel;
        % all parameters are used in each model
        myForwardModels(i).PMap  = 1:ndim; 
        
    end

    
    %% Likelihood description
    % Here, the discrepancy for each data structure |y| are chosen to be
    % independent and identically distributed Gaussian random variables.
    % For the current case, 2*standard deviations of the experimental
    % measurements is chosen as the prior.

    DiscrepancyOpts(i).Type = 'Gaussian';
    DiscrepancyOpts(i).Parameters = 1e-1;    
    
    % DiscrepancyPriorOpts1.Name = 'Prior of sigma 1';
    % DiscrepancyPriorOpts1.Marginals(1).Name = 'Sigma1';
    % DiscrepancyPriorOpts1.Marginals(1).Type = 'Uniform';
    % DiscrepancyPriorOpts1.Marginals(1).Parameters = [0, 2*std(output_raw.Fy03)];
    % DiscrepancyPrior1 = uq_createInput(DiscrepancyPriorOpts1);
    %
    % DiscrepancyOpts(1).Type = 'Gaussian';
    % DiscrepancyOpts(1).Prior = DiscrepancyPrior1;
    
    %
    % DiscrepancyPriorOpts2.Name = 'Prior of sigma 2';
    % DiscrepancyPriorOpts2.Marginals(1).Name = 'Sigma2';
    % DiscrepancyPriorOpts2.Marginals(1).Type = 'Uniform';
    % DiscrepancyPriorOpts2.Marginals(1).Parameters = [0, 2*std(output_raw.Fy05)];
    % DiscrepancyPrior2 = uq_createInput(DiscrepancyPriorOpts2);
    %
    % DiscrepancyOpts(2).Type = 'Gaussian';
    % DiscrepancyOpts(2).Prior = DiscrepancyPrior2;
    %
    % DiscrepancyPriorOpts3.Name = 'Prior of sigma 3';
    % DiscrepancyPriorOpts3.Marginals(1).Name = 'Sigma3';
    % DiscrepancyPriorOpts3.Marginals(1).Type = 'Uniform';
    % DiscrepancyPriorOpts3.Marginals(1).Parameters = [0, 2*std(output_raw.Fy08)];
    % DiscrepancyPrior3 = uq_createInput(DiscrepancyPriorOpts3);
    %
    % DiscrepancyOpts(3).Type = 'Gaussian';
    % DiscrepancyOpts(3).Prior = DiscrepancyPrior3;
    %
    % DiscrepancyPriorOpts4.Name = 'Prior of sigma 4';
    % DiscrepancyPriorOpts4.Marginals(1).Name = 'Sigma4';
    % DiscrepancyPriorOpts4.Marginals(1).Type = 'Uniform';
    % DiscrepancyPriorOpts4.Marginals(1).Parameters = [0, 2*std(output_raw.Fy10)];
    % DiscrepancyPrior4 = uq_createInput(DiscrepancyPriorOpts4);
    %
    % DiscrepancyOpts(4).Type = 'Gaussian';
    % DiscrepancyOpts(4).Prior = DiscrepancyPrior4;
    
    % DiscrepancyOpts.Type = 'Gaussian';
    % DiscrepancyOpts.Parameters = 1e-6;
    
    
end




%% check location of ECNAeroModule
% path_found  = findAeroModulePath();

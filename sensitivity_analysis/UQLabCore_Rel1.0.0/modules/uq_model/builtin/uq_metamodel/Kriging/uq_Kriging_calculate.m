function success = uq_Kriging_calculate(current_model)
% UQ_KRIGING_CALCULATE(CURRENT_MODEL): calculate
% the Kriging metamodel specified in CURRENT_MODEL
%
% See also: UQ_KRIGING_OPTIMIZER, UQ_KRIGING_EVAL_F

success = 0; 

%% argument and consistency checks
% let's check the model is of type "uq_metamodel"
if ~strcmp(current_model.Type, 'uq_metamodel')
    error('Error: uq_metamodel cannot handle objects of type %s', current_model.Type);
end

%% Reporting
DisplayLevel = current_model.Internal.Display ;
if DisplayLevel
    fprintf('---   Calculating the Kriging metamodel...                              ---\n')
end

%% Generate the initial experimental design
% Get X
[current_model.ExpDesign.X, current_model.ExpDesign.U] = uq_getExpDesignSample(current_model);
% Get Y
current_model.ExpDesign.Y = uq_eval_ExpDesign(current_model,current_model.ExpDesign.X);
% Update the number of output variables of the model and store it
Nout = size(current_model.ExpDesign.Y, 2);
current_model.Internal.Runtime.Nout = Nout;


%% ----  Kriging entry point
X = current_model.ExpDesign.U;
try
    if isfield(current_model.Internal,'Input') && ...
            (isfield(current_model.Internal.Input,'nonConst') || ...
            isprop(current_model.Internal.Input,'nonConst'))
        nonConst = current_model.Internal.Input.nonConst;
    end
    if isfield(current_model.Internal,'nonConst')
        nonConst = current_model.Internal.nonConst;
    end
    X = X(:,nonConst);
catch
    nonConst = 1:size(X,2);
end
N = size(X,1) ;
current_model.Internal.Runtime.N = N;
%% calculate F and store it
current_model.Internal.Runtime.current_output = 1;
evalF_handle = current_model.Internal.Kriging(1).Trend.Handle ;
F = evalF_handle( X, current_model ) ;
current_model.Internal.Kriging.Trend.F = F;

%% get the estimation method 
estim_method = current_model.Internal.Kriging.GP.EstimMethod;

%% parse the optimization options
optim_options = uq_Kriging_initialize_optimizer(current_model);
% Cycle through each output
for oo = 1 : Nout
    % copy the necessary information about the Kriging options to the various
    % output coordinates 
    if oo > 1
        current_model.Internal.Kriging(oo) =  current_model.Internal.Kriging(1);
    end
    
    % store the current output
    current_model.Internal.Runtime.current_output = oo;
    
    % Get the current Y
    Y = current_model.ExpDesign.Y(:,oo) ;
    
    %% find optimal theta
    [theta,Jstar, fcount, nIter, exitflag] = ...
        uq_Kriging_optimizer(X, Y, optim_options, current_model) ;
    % store optimization-related results inside current_model
    % optimal theta
    current_model.Internal.Kriging(oo).Optim.Theta = theta; 
    % objective function value at optimal theta
    current_model.Internal.Kriging(oo).Optim.ObjFun = Jstar; 
    % number of function evaluations
    current_model.Internal.Kriging(oo).Optim.nEval = fcount; 
    % number of iterations
    current_model.Internal.Kriging(oo).Optim.nIter = nIter ; 
    % exit flag that describes the exit condition of the optimization process (see
    % uq_Kriging_optimizer for more details)
    current_model.Internal.Kriging(oo).Optim.ExitFlag = exitflag ;
    
    %% calculate R 
    evalR_handle = current_model.Internal.Kriging(oo).GP.Corr.Handle ;
    R = evalR_handle( X, X, theta, current_model.Internal.Kriging(oo).GP.Corr);
    %  and store it
    current_model.Internal.Kriging(oo).GP.R = R ;
    
    %% calculate sigma^2 and beta 
    if strcmpi(estim_method,'cv')
        CV_K = current_model.Internal.Kriging(oo).GP.CV.CV_K;
        [ errors,  additional_metrics, auxMatrices] = ...
        uq_Kriging_calc_leaveKout(CV_K, N , Y, R, F, 'default' ) ;
        parameters.errors = errors;
        parameters.sigma_i_Sq = additional_metrics.sigma_i_Sq;
        parameters.CV_K = CV_K;
           % calculate sigma^2 and error metrics and store them
        sigmaSQ = uq_Kriging_calc_sigmaSq( parameters, estim_method );                    
        current_model.Internal.Kriging(oo).GP.sigmaSQ = sigmaSQ; 
                    
        % calculate beta and store it
        beta = uq_Kriging_calc_beta( ...
                        F, ...
                        current_model.Internal.Kriging(oo).Trend.Type, ...
                        Y, ...
                        'standard', ...
                        auxMatrices );
        current_model.Internal.Kriging(oo).Trend.beta = beta;
    else % currently only 'ml' is a viable alternative
       parameters.Y = Y;
       parameters.N = N;
       parameters.F = F;
       
       [ errors,  additional_metrics, auxMatrices] = ...
            uq_Kriging_calc_leaveKout(N, N , Y, R, F, 'ml_estimation' ) ;
        
        % calculate beta and store it
        beta = uq_Kriging_calc_beta( ...
                        F, ...
                        current_model.Internal.Kriging(oo).Trend.Type, ...
                        Y, ...
                        'qr', ...
                        auxMatrices );
        current_model.Internal.Kriging(oo).Trend.beta = beta;
        
        % calculate sigma^2 and error metrics and store them
        parameters.beta = beta;
        if ~isnan(auxMatrices.cholR)
            parameters.Ytilde = auxMatrices.Ytilde;
            parameters.Ftilde = auxMatrices.Ftilde;
            sigmaSQ = uq_Kriging_calc_sigmaSq( parameters, 'ml_chol' );
        else
            parameters.Rinv = auxMatrices.Rinv;
            sigmaSQ = uq_Kriging_calc_sigmaSq( parameters, 'ml_NOchol' );
        end
        current_model.Internal.Kriging(oo).GP.sigmaSQ = sigmaSQ; 
    end
    %% store error metrics
    varY = var(Y, 1, 1);                     
    current_model.Error(oo).LOO =  mean(cell2mat(errors))/varY;    
    current_model.Internal.Error(oo).varY = varY ;
    current_model.Internal.Error(oo).LOOmean = additional_metrics.LOOmean;
    current_model.Internal.Error(oo).LOOsd = additional_metrics.LOOsd;
    current_model.Internal.Error(oo).sigma_i_Sq = additional_metrics.sigma_i_Sq;

    %% store cache if requested to do so
    if current_model.Internal.KeepCache
        % store all the auxiliary matrices that can be useful for speeding
        % up predictions
        current_model.Internal.Kriging(oo).Cached.cholR  = auxMatrices.cholR;
        current_model.Internal.Kriging(oo).Cached.Rinv   = auxMatrices.Rinv;
        current_model.Internal.Kriging(oo).Cached.FTRinv = auxMatrices.FTRinv;
        current_model.Internal.Kriging(oo).Cached.FTRinvF= auxMatrices.FTRinvF;
    else
        current_model.Internal.Kriging(oo).Cached = [] ;
    end

    %% store a copy of the most important Kriging results inside current_model.Kriging
    current_model.Kriging(oo).beta = current_model.Internal.Kriging(oo).Trend.beta;
    current_model.Kriging(oo).sigmaSQ = current_model.Internal.Kriging(oo).GP.sigmaSQ;
    current_model.Kriging(oo).theta = current_model.Internal.Kriging(oo).Optim.Theta;
end % end of kriging per output loop


% Raise the flag that the metamodel  has been calculated
current_model.Internal.Runtime.isCalculated = 1;

if DisplayLevel
    fprintf('---   Calculation finished!                                             ---\n')
end

success = 1;

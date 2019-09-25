function [theta,Jstar, fcount, nIter, exitflag] = ...
    uq_Kriging_optimizer( X, Y, optim_options, current_model )
% [theta,Jstar, fcount, nIter, exitflag]=UQ_KRIGING_OPTIMIZER( X, Y, optim_options, ...
% current_model) : returns the optimal values of the hyperparameters (theta) together 
% with some information about the optimization that took place:
%   Jstar: the value of the objective function on the optimum that was
%   found
%   fcount: number of objective function evaluations that took place
%   nIter: the number of iterations/generations of the optimization
%   algorithm
%   exitflag: signals the reason why the optimization stopped. The values
%   of the exitflag are depending on the optmization method that is used
%   (see e.g. fmincon, ga, etc.)
%
% See also UQ_KRIGING_CALCULATE, UQ_KRIGING_INITIALIZE_OPTIMIZER, 
% UQ_KRIGING_EVAL_J_OF_THETA_CV, UQ_KRIGING_EVAL_J_OF_THETA_ML


%% READ options from the current module
% Obtain the current output
current_output = current_model.Internal.Runtime.current_output ;
% Retrieve the handle of the appropriate objective function
switch lower(current_model.Internal.Kriging(current_output).GP.EstimMethod)
    case 'ml'
        objFunHandle = str2func('uq_Kriging_eval_J_of_theta_ML') ;
    case 'rml'
        objFunHandle = str2func('uq_Kriging_eval_J_of_theta_RML') ;
    case 'cv'
        objFunHandle = str2func('uq_Kriging_eval_J_of_theta_CV') ;
end

M = current_model.Internal.Runtime.M ;

%% Organize in a simple structure the necessary parameters in order to solve
%  this specific optimization problem with these specific options
parameters.X = X;
parameters.Y = Y;
parameters.N = current_model.Internal.Runtime.N ;
parameters.F = current_model.Internal.Kriging(current_output).Trend.F;
parameters.trend_type = current_model.Internal.Kriging(current_output).Trend.Type;
parameters.CorrOptions = current_model.Internal.Kriging(current_output).GP.Corr ;

estim_method = current_model.Internal.Kriging(current_output).GP.EstimMethod;
switch lower(estim_method)
    case 'cv'
        parameters.CV_K = current_model.Internal.Kriging(current_output).GP.CV.CV_K;
    case 'ml'
end
%% Find the size of the optimization variable

% First lets see whether the build-in uq_Kriging_eval_R is used to calculate R:
if strcmpi(func2str(current_model.Internal.Kriging(current_output).GP.Corr.Handle),...
        'uq_eval_kernel')
    % in this case the build-in correlation matrix calculator is used
    % Keeping track of the non-constants:
    nonConstIdx = current_model.Internal.Runtime.nonConstIdx;
    % isotropy of the correlation function
    isIsotropic = current_model.Internal.Kriging(current_output).GP.Corr.Isotropic;
    
    % Determine the number of optimization variables
    % NOTE: It is assumed for anisotropic correlation functions there exists
    % *one* optimization variable per dimension. That might not be the case in
    % a future release.
    if isIsotropic
        % isotropic correlation function
        nVars = 1;
        % do nothing about constants (assume that at least one variable is
        % non-constant)
    else
        nVars = length(nonConstIdx);
    end
else
    % Keeping track of the non-constants:
    nonConstIdx = current_model.Internal.Runtime.nonConstIdx;
    
    % in this case a user-defined correlation matrix calculator is used
    % lets try to find the optimization variable size by the initial value
    % and/or the bounds that were defined (if both are defined they should be 
    % consistent!)
    
    nVars = [];
    if isfield(current_model.Internal.Kriging(current_output).Optim, 'InitialValue')
        initval_definition = current_model.Internal.Kriging(current_output).Optim.InitialValue;
        nVars = [nVars; length(initval_definition)];
    end
   
    if isfield(current_model.Internal.Kriging(current_output).Optim, 'Bounds')
        % Keep only  the bounds that correspond to the varying dimensions
        % of the metamodel
        bounds_definition = current_model.Internal.Kriging(current_output).Optim.Bounds;
        nVars = [nVars; size(bounds_definition,2)];
    end
    
    switch length(nVars)
        case 0 
            error('Kriging optimization error: either the InitialValue or Bounds need to be defined!')
        case 1
            % do nothing
        case 2
            if nVars(1)~=nVars(2)
                error('Kriging optimization error: Inconsistent dimensions between the InitialValue and Bounds!')
            end
            nVars = nVars(1);
    end
    
end


%% Execute optimization based on the selected method
switch lower(current_model.Internal.Kriging(current_output).Optim.Method)
    case 'none'
        theta = current_model.Internal.Kriging(current_output).Optim.InitialValue;
        Jhandle = @(theta)objFunHandle(theta, parameters) ;
        Jstar = Jhandle(theta);
        fcount = 1;
        nIter = 0;
        exitflag = 1;
    case {'gradbased','lbfgs', 'bfgs'}
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        theta0 = current_model.Internal.Kriging(current_output).Optim.InitialValue;
        current_model.Internal.Kriging(current_output).Optim.InitialObjFun = ...
            objFunHandle(theta0, parameters); %initial value of objective function
        %Gradient-based optimization : currently using fmincon
        
        [theta, Jstar, exitflag, output] = ...
            fmincon(@(theta)objFunHandle(theta, parameters), ...
            theta0,...
            [], [], [], [], LB, UB,[], optim_options) ;
        
        fcount = output.funcCount ;
        nIter = output.iterations ;
        
    case 'knitro'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        theta0 = current_model.Internal.Kriging(current_output).Optim.InitialValue;
        current_model.Internal.Kriging(current_output).Optim.InitialObjFun = ...
            objFunHandle(theta0, X, current_model); %initial value of objective function

        
        [theta, Jstar, exitflag, output] = ...
            ktrlink(@(theta)objFunHandle(theta, parameters),...
            theta0,[],[],[],[], LB, UB,[],optim_options) ;
        fcount = output.funcCount ;
        nIter = output.iterations ;
    case 'ga'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        % Vanilla Genetic Algorithm optimization

        [theta, Jstar, exitflag, output] = ...
            ga(@(theta)objFunHandle(theta, parameters), ...
            nVars,[], [], [], [], LB, UB, [], optim_options ) ;
        fcount = output.funccount ;
        nIter = output.generations ;
    case 'hga'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        % Hybrid Genetic Algorithm optimization
        
        
        [theta, ~, exitflag.GA, output] = ga(@(theta)objFunHandle(theta, parameters), ...
            nVars,...
            [], [], [], [], LB, UB, [], optim_options.ga ) ;
        fcountGA = output.funccount ;
        nIterGA = output.generations ;
        % Refine the result by executing a gradient based optimization
        [theta, Jstar, exitflag.BFGS, output] = fmincon(@(theta)objFunHandle(theta, parameters), ...
            theta,...
            [], [], [], [], LB, UB,[], optim_options.grad) ;
        
        fcountGRAD = output.funcCount ;
        nIterGRAD = output.iterations ;
        
        fcount = fcountGA + fcountGRAD ;
        %nIterGRAD is not taken into account when calculating the total
        %number of iterations 
        nIter = nIterGA; 
        
    case 'sade'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        % Self-Adaptive Differential Evolution
        Npop = current_model.Internal.Kriging(current_output).Optim.SADE.nPop ;
        
        [ theta, Jstar, output ] = ...
            uq_optim_sade(@(theta)objFunHandle(theta, parameters)...
            , nVars, Npop, LB, UB, optim_options);
        
        fcount = output.fcount ;
        nIter = output.niter ;
        exitflag = output.exitflag ;
        % store some extra interesting fields
        current_model.Internal.Kriging(current_output).Optim.Strategies = output.strategies ;
        current_model.Internal.Kriging(current_output).Optim.pStrategies = output.pStrategies ;
        current_model.Internal.Kriging(current_output).Optim.CRm = output.CRm ;
    case 'hsade'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        Npop = current_model.Internal.Kriging(current_output).Optim.HSADE.nPop ;
        
        % Differential Evolution
        [ theta, Jstar, output ] = ...
            uq_optim_sade(@(theta)objFunHandle(theta, parameters)...
            , nVars, Npop, LB, UB, optim_options.sade);
        fcountDE = output.fcount ;
        nIterDE = output.niter ;
        exitflag.SADE = output.exitflag ;
        
        % store some extra interesting fields
        current_model.Internal.Kriging(current_output).Optim.Strategies = output.strategies ;
        current_model.Internal.Kriging(current_output).Optim.pStrategies = output.pStrategies ;
        current_model.Internal.Kriging(current_output).Optim.CRm = output.CRm ;
        
        % Refine the result my executing a gradient based optimization
        [theta, Jstar, exitflag.BFGS, output] = fmincon(@(theta)objFunHandle(theta, parameters), ...
            theta,...
            [], [], [], [], LB, UB,[], optim_options.grad) ;
        
        fcountGRAD = output.funcCount ;
        nIterGRAD = output.iterations ;
        
        fcount = fcountDE + fcountGRAD ;
        %nIterGRAD is not taken into account when calculating the total
        %number of iterations
        nIter = nIterDE;
 
            case 'cmaes'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        % Vanilla Genetic Algorithm optimization

        [theta, Jstar, exitflag, output] = ...
            uq_cmaes(@(theta)objFunHandle(theta, parameters), ...
            [],[], LB, UB, optim_options ) ;
        fcount = output.funccount ;
        nIter = output.iterations ;
        
    case 'hcmaes'
        LB = current_model.Internal.Kriging(current_output).Optim.Bounds(1,:);
        UB = current_model.Internal.Kriging(current_output).Optim.Bounds(2,:);
        % Hybrid CMA-ES optimization
        
        [theta, Jstar, exitflag.CMAES, output] = ...
            uq_cmaes(@(theta)objFunHandle(theta, parameters), ...
            [],[], LB, UB, optim_options.cmaes ) ;        

        fcountCMAES = output.funccount ;
        nIterCMAES = output.iterations ;
        
        % Refine the result by executing a gradient based optimization
        [theta, Jstar, exitflag.BFGS, output] = fmincon(@(theta)objFunHandle(theta, parameters), ...
            theta,...
            [], [], [], [], LB, UB,[], optim_options.grad) ;
        
        fcountGRAD = output.funcCount ;
        nIterGRAD = output.iterations ;
        
        fcount = fcountCMAES + fcountGRAD ;
        %nIterGRAD is not taken into account when calculating the total
        %number of iterations 
        nIter = nIterCMAES; 
        
    otherwise
        error('Error: Unknown method for calculating Kriging covariance matrix hyperparameters!')
end





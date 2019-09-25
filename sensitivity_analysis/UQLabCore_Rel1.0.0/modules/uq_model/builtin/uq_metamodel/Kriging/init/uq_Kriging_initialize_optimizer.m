function optim_options = uq_Kriging_initialize_optimizer(current_model)
% optim_options = UQ_KRIGING_INITIALIZE_OPTIMIZER(KRGModel) Parses various optimization 
% options depending on the optimization method that is selected. It is assumed that 
% for each output the same optimization method will be used
%
% See also: UQ_KRIGING_OPTIMIZER, UQ_KRIGING_CALCULATE


nonConstIdx = current_model.Internal.Runtime.nonConstIdx;

%% 
if strcmpi(func2str(current_model.Internal.Kriging(1).GP.Corr.Handle),...
        'uq_kriging_eval_r')
    
    %% filter out constants from the initial value and bounds definition
    isIsotropic = current_model.Internal.Kriging(1).GP.Corr.Isotropic;
    if isIsotropic
        % do nothing about constants (assume that at least one variable is
        % non-constant)
    else
        % anisotropic correlation function
        if isfield(current_model.Internal.Kriging(1).Optim, 'InitialValue')
            current_model.Internal.Kriging(1).Optim.InitialValue = ...
             current_model.Internal.Kriging(1).Optim.InitialValue(nonConstIdx);
        end

        if isfield(current_model.Internal.Kriging(1).Optim, 'Bounds')
            % Keep only  the bounds that correspond to the varying dimensions
            % of the metamodel
            current_model.Internal.Kriging(1).Optim.Bounds = ...
                current_model.Internal.Kriging(1).Optim.Bounds(:,nonConstIdx);
        end
    end
else
   isIsotropic = false; 
end

%% parse optimization options
switch lower(current_model.Internal.Kriging(1).Optim.Method)
    case 'none'
       optim_options = []; 
       return
    case {'gradbased','lbfgs', 'bfgs'}
       optim_options = optimset(...
    'Display', lower(current_model.Internal.Kriging(1).Optim.Display),...
    'MaxIter', current_model.Internal.Kriging(1).Optim.MaxIter,...
    'Algorithm','interior-point',...
    'Hessian',{'lbfgs',current_model.Internal.Kriging(1).Optim.BFGS.nLM},...
    'AlwaysHonorConstraints','none',...
    'TolFun', current_model.Internal.Kriging(1).Optim.Tol );
    case 'ga'    
        optim_options = gaoptimset(...
        'Display',lower(current_model.Internal.Kriging(1).Optim.Display), ...
        'Generations',current_model.Internal.Kriging(1).Optim.MaxIter,...
        'PopulationSize', current_model.Internal.Kriging(1).Optim.GA.nPop,...
        'StallGenLimit', current_model.Internal.Kriging(1).Optim.GA.nStall,...
        'TolFun', current_model.Internal.Kriging(1).Optim.Tol);
    case 'hga'
        optim_options.ga = ...
            gaoptimset('Display',lower(current_model.Internal.Kriging(1).Optim.Display), ...
            'Generations',current_model.Internal.Kriging(1).Optim.MaxIter,...
            'PopulationSize', current_model.Internal.Kriging(1).Optim.HGA.nPop,...
            'StallGenLimit', current_model.Internal.Kriging(1).Optim.HGA.nStall,...
            'TolFun', current_model.Internal.Kriging(1).Optim.Tol);
        optim_options.grad = ...
            optimset('Display',lower(current_model.Internal.Kriging(1).Optim.Display),...
            'MaxIter', current_model.Internal.Kriging(1).Optim.MaxIter,...
            'Algorithm','interior-point', 'Hessian',{'lbfgs',current_model.Internal.Kriging(1).Optim.HGA.nLM},...
            'AlwaysHonorConstraints','none',...
            'TolFun', current_model.Internal.Kriging(1).Optim.Tol );
        
    case 'sade'
        optim_options.Display = lower(current_model.Internal.Kriging(1).Optim.Display) ;
        optim_options.MaxIter = current_model.Internal.Kriging(1).Optim.MaxIter ;
        optim_options.TolFun = current_model.Internal.Kriging(1).Optim.Tol;
        optim_options.nStall = current_model.Internal.Kriging(1).Optim.SADE.nStall ;
        optim_options.Strategies = current_model.Internal.Kriging(1).Optim.SADE.Strategies ;
        optim_options.pStr = current_model.Internal.Kriging(1).Optim.SADE.pStr ;
        optim_options.CRm = current_model.Internal.Kriging(1).Optim.SADE.CRm ;
    case 'hsade'
        optim_options.sade.Display = lower(current_model.Internal.Kriging(1).Optim.Display) ;
        optim_options.sade.MaxIter = current_model.Internal.Kriging(1).Optim.MaxIter ;
        optim_options.sade.TolFun = current_model.Internal.Kriging(1).Optim.Tol;
        optim_options.sade.nStall = current_model.Internal.Kriging(1).Optim.HSADE.nStall ;
        optim_options.sade.Strategies = current_model.Internal.Kriging(1).Optim.HSADE.Strategies ;
        optim_options.sade.pStr = current_model.Internal.Kriging(1).Optim.HSADE.pStr ;
        optim_options.sade.CRm = current_model.Internal.Kriging(1).Optim.HSADE.CRm ;
        
        optim_options.grad = ...
            optimset('Display',lower(current_model.Internal.Kriging(1).Optim.Display),...
            'MaxIter', current_model.Internal.Kriging(1).Optim.MaxIter,...
            'Algorithm','interior-point',...
            'Hessian',{'lbfgs',current_model.Internal.Kriging(1).Optim.HSADE.nLM},...
            'AlwaysHonorConstraints','none',...
            'TolFun', current_model.Internal.Kriging(1).Optim.Tol );
        
    case 'cmaes'
        optim_options.Display = lower(current_model.Internal.Kriging(1).Optim.Display) ;
        optim_options.MaxIter = current_model.Internal.Kriging(1).Optim.MaxIter ;
        optim_options.TolX = current_model.Internal.Kriging(1).Optim.Tol ;

        optim_options.lambda = current_model.Internal.Kriging(1).Optim.CMAES.nPop ;
        optim_options.nStallMax = current_model.Internal.Kriging(1).Optim.CMAES.nStall ;
        
        % Output is not vectorized:
        optim_options.isVectorized = false ;
        
    case 'hcmaes'

        optim_options.cmaes.Display = lower(current_model.Internal.Kriging(1).Optim.Display) ;
        optim_options.cmaes.MaxIter = current_model.Internal.Kriging(1).Optim.MaxIter ;
        optim_options.cmaes.TolX = current_model.Internal.Kriging(1).Optim.Tol ;
        
        optim_options.cmaes.lambda = current_model.Internal.Kriging(1).Optim.HCMAES.nPop ;
        optim_options.cmaes.nStallMax = current_model.Internal.Kriging(1).Optim.HCMAES.nStall ;
        
        % Output is not vectorized:
        optim_options.cmaes.isVectorized = false ;
        
        optim_options.grad = ...
            optimset('Display',lower(current_model.Internal.Kriging(1).Optim.Display),...
            'MaxIter', current_model.Internal.Kriging(1).Optim.MaxIter,...
            'Algorithm','interior-point', 'Hessian',{'lbfgs',current_model.Internal.Kriging(1).Optim.HCMAES.nLM},...
            'AlwaysHonorConstraints','none',...
            'TolFun', current_model.Internal.Kriging(1).Optim.Tol );
end


function success = uq_Kriging_initialize( current_model )
% SUCCESS = UQ_KRIGING_INITIALIZE(KRGMODEL): Initialize a Kriging model based on the
%     user-specified options.
% 
% See also: UQ_KRIGING_INITIALIZE_CUSTOM, UQ_PCE_INITIALIZE,
% UQ_INITIALIZE_UQ_METAMODEL

success = 0;

M = current_model.Internal.Runtime.M;

%% UNPROCESSED FIELDS
skipFields = {'Type','Name','MetaType','Input','FullModel','ExpDesign', 'Display','ValidationSet'};

%% DEFAULT VALUES
% ExpDesign
ScalingDefaults = 1 ;
keepCache = 1; %flag to keep the cached matrices after the metamodel has been calculated

% Trend
TrendDefaults.Type = 'ordinary' ;
TrendDefaults.Degree = 0 ;
TrendDefaults.TruncOptions.qNorm = 1;
TrendDefaults.CustomF = [];
TrendDefaults.PolyTypes = cell(1,M);
[TrendDefaults.PolyTypes{1:M}] = deal('simple_poly' );
TrendDefaults.Handle = @uq_Kriging_eval_F;

% Gaussian Process
CorrDefaults.Family = 'matern-5_2';
CorrDefaults.Type = 'ellipsoidal';
CorrDefaults.Isotropic = false;
CorrDefaults.Nugget = 1e-10;
CorrDefaults.Handle = @uq_eval_Kernel;

EstMethodDefaults.EstimMethod = 'CV' ;
EstMethodDefaults.CV.LeaveKOut = 1;

% Optimization
OptimFieldsDatatypes.INITIALVALUE = 'double';
OptimFieldsDatatypes.BOUNDS = 'double';
OptimFieldsDatatypes.MAXITER = 'double';
OptimFieldsDatatypes.TOL = 'double';
OptimFieldsDatatypes.DISPLAY = 'char';
OptimFieldsDatatypes.NLM = 'double';
OptimFieldsDatatypes.NPOP = 'double';
OptimFieldsDatatypes.NSTALL = 'double';
OptimFieldsDatatypes.STRATEGIES = 'cell';
OptimFieldsDatatypes.PSTR = 'double' ;
OptimFieldsDatatypes.CRM = 'double';

% Known optimization methods
OptimKnownMethods = {'None','GA', 'LBFGS','BFGS','HGA', 'SADE', 'HSADE',...
    'CMAES','HCMAES','CE','HCE','1p1CMAES','H1p1CMAES'} ;
% Known optimization methods with respect to toolbox
OptimKnownOptimToolboxMethods = {'LBFGS','BFGS','HSADE'} ;
OptimKnownGlobalOptimToolboxMethods = 'GA' ;
OptimKnownOptimAndGlobalToolboxMethods = 'HGA' ;
% Required field for each method:
% NOTE: The required fields for the SPECIFIC method are set below in OptimDefaults_Method
OptimReqFields.NONE = {'InitialValue'};
OptimReqFields.GA   = {'Bounds', 'MaxIter', 'Tol', 'Display'} ;
OptimReqFields.SADE = OptimReqFields.GA ; % Subject to change in a future release
OptimReqFields.BFGS   = {'InitialValue','Bounds', 'MaxIter', 'Tol', 'Display'} ;
OptimReqFields.HGA = OptimReqFields.GA;
OptimReqFields.HSADE = OptimReqFields.SADE;
OptimReqFields.KNITRO = OptimReqFields.BFGS ; % Subject to change in a future release
OptimReqFields.CMAES   = {'Bounds', 'MaxIter', 'Tol', 'Display'} ;
OptimReqFields.CE   = {'Bounds', 'MaxIter', 'Tol', 'Display'} ;
OptimReqFields.HCMAES = OptimReqFields.CMAES;

% Method-specific default values
OptimDefaults_Method.NONE = [];
OptimDefaults_Method.BFGS.nLM = 5 ;
OptimDefaults_Method.GA.nPop =  30; 
OptimDefaults_Method.GA.nStall = 5;
OptimDefaults_Method.SADE = OptimDefaults_Method.GA ;% Subject to change in a future release
OptimDefaults_Method.SADE.Strategies = {'rand_1_bin', 'rand_2_bin', ...
    'rand_to_best_2_bin','curr_to_rand_1',...
    'best_1_bin','rand_2_bin','rand_to_best_2_bin'} ;
% keep the number of strategies
nStrategies = length(OptimDefaults_Method.SADE.Strategies) ;
OptimDefaults_Method.SADE.pStr = (1/nStrategies)* ones(nStrategies,1) ;
OptimDefaults_Method.SADE.CRm = repmat(0.5,nStrategies,1);
OptimDefaults_Method.HGA = merge_structures(OptimDefaults_Method.BFGS,...
    OptimDefaults_Method.GA) ;
OptimDefaults_Method.HSADE = merge_structures(OptimDefaults_Method.BFGS,...
    OptimDefaults_Method.SADE) ;
OptimDefaults_Method.KNITRO = OptimDefaults_Method.BFGS ;

OptimDefaults_Method.CMAES.nPop =  30;
OptimDefaults_Method.CMAES.nStall = 5;

OptimDefaults_Method.HCMAES = merge_structures(OptimDefaults_Method.BFGS,...
    OptimDefaults_Method.CMAES) ;
% Default values when nothing is set by the user
OptimDefaults.InitialValue = 1 ;
OptimDefaults.Bounds = [ 1e-3 ; 10 ]; % [LOWER ; UPPER]
OptimDefaults.Method = 'HGA' ;
OptimDefaults.HGA = OptimDefaults_Method.HGA ;
OptimDefaults.MaxIter = 20 ;
OptimDefaults.Tol = 1e-4 ;
OptimDefaults.Display = 'none' ;

%% RETRIEVE THE OPTIONS AND PARSE THEM
Options = current_model.Options;

%% PARSE GLOBAL DISPLAY LEVEL
% Get the global verbosity level
DisplayLevel = current_model.Internal.Display;

% If Optim.Display is not manually set update it based on the global
% verbosity level
OptimDisplay_EXISTS = isfield(Options, 'Optim') && isfield(Options.Optim, 'Display') &&...
    ~isempty(Options.Optim.Display);
OptimMethod_NONE = isfield(Options, 'Optim') && isfield(Options.Optim, 'Method') &&...
    ~isempty(Options.Optim.Method) && ...
    strcmpi(Options.Optim.Method, 'none');

switch lower(DisplayLevel)
    case 0
        if ~OptimDisplay_EXISTS && ~OptimMethod_NONE
           Options.Optim.Display = 'none'; 
        end
    case 1
        if ~OptimDisplay_EXISTS && ~OptimMethod_NONE
           Options.Optim.Display = 'final'; 
        end
    case 2
        if ~OptimDisplay_EXISTS && ~OptimMethod_NONE
           Options.Optim.Display = 'iter'; 
        end
    otherwise
        EVT.Type = 'W';
        EVT.Message = sprintf('Unknown display option: %s. Using the default value instead.', ...
            num2str(DisplayLevel));
        EVT.eventID = 'uqlab:metamodel:kriging:init:display_invalid';
        uq_logEvent(current_model, EVT);
        % Set the default display level
        DisplayLevel = 1;
end

%% METATYPE
% Meta Type
if ~isfield(Options, 'MetaType') || isempty(Options.MetaType)
    error('MetaType must be specified.');
end

uq_addprop(current_model, 'MetaType', Options.MetaType);
uq_addprop(current_model, 'Internal');

%% INPUT 
% Check whether an INPUT object has been defined
INPUT_EXISTS = isfield(current_model.Internal,'Input') && ...
    ~isempty(current_model.Internal.Input);
% Temporary fix for the case that both an ED and Input have been specified: 
if ~INPUT_EXISTS && isfield(Options,'Input') && ~isempty(Options.Input)
   INPUT_EXISTS = true;
   current_model.Internal.Input = Options.Input;
end

% When an input module is specified it can be due to one (or more) of the
% following reasons:
% - An experimental design needs to be generated according to the
% probability distribution of the INPUT
% - A special type of scaling needs to take place (by isoprobabilistic
% transform)
if INPUT_EXISTS
    if any(strcmpi(current_model.ExpDesign.Sampling, {'user', 'data'}))
        % do a consistency check of the input dimension
        if length(current_model.Internal.Input.Marginals) ~= ...
                size(current_model.ExpDesign.X, 2)
            error('Input dimension inconsistency!')
        end
    end
end

%% SCALING (Auxiliary space)
% The struct case is included for the case of having a user defined
% auxiliary space (not yet implemented!)
[scale, Options] = uq_process_option(Options, 'Scaling',...
    ScalingDefaults, {'struct','logical','double', 'uq_input'});
if scale.Invalid
    EVT.Type = 'W';
    EVT.Message = 'The Scaling option was invalid. Using the default value instead.';
    EVT.eventID = 'uqlab:metamodel:kriging:init:scaling_invalid';
    uq_logEvent(current_model, EVT);
end
if scale.Missing
    %do nothing, just silently assign the default value
    EVT.Type = 'D';
    EVT.Message = 'The Scaling option was missing. Assigning the default value.';
    EVT.eventID = 'uqlab:metamodel:kriging:init:scaling_missing';
    uq_logEvent(current_model, EVT);
end

% get and store the Scaling value
current_model.Internal.Scaling = scale.Value;
SCALING = scale.Value ;
SCALING_BOOL = isa(SCALING, 'double') || isa(SCALING, 'logical') || isa(SCALING, 'int');

if SCALING_BOOL && SCALING
    if INPUT_EXISTS
        % scale the data as U = (X - muX)/stdX where muX,stdX are computed from the specified input distribution
        input_moments = reshape([current_model.Internal.Input.Marginals(:).Moments],2,[]);
        current_model.Internal.ExpDesign.muX = input_moments(1,:); % mean
        current_model.Internal.ExpDesign.sigmaX = input_moments(2,:);% standard deviation
    else
        % scale the data as U = (X - muX)/stdX where muX,stdX are computed from the available data
        current_model.Internal.ExpDesign.muX = mean(current_model.ExpDesign.X);
        current_model.Internal.ExpDesign.sigmaX = std(current_model.ExpDesign.X);
    end
end

if ~SCALING_BOOL
    if ~INPUT_EXISTS
       error('An Input object needs to be specified for the scaling option selected!') 
    end
    if isa(SCALING, 'uq_input')
        % do nothing
    else
       current_model.Internal.Scaling = uq_createInput(SCALING, '-private');
    end
end


%% Keep cache?
[ckeep, Options] = uq_process_option(Options, 'KeepCache',...
    keepCache, {'logical','double'});
if ckeep.Invalid 
    EVT.Type = 'W';
    EVT.Message = 'The KeepCache option was invalid. Using the default value instead.';
    EVT.eventID = 'uqlab:metamodel:kriging:init:keepcache_invalid';
    uq_logEvent(current_model, EVT);
end
if ckeep.Missing
    %do nothing, just silently assign the default value
end
current_model.Internal.KeepCache = ckeep.Value ;

%% Parse Trend-Related Options
if isfield(Options, 'Trend') && ~isempty(Options.Trend) && ...
    isfield(Options.Trend, 'Handle') && ~isempty(Options.Trend.Handle)
    trend_handle = Options.Trend.Handle ;
else 
    trend_handle = TrendDefaults.Handle ;
end

if strcmp(char(trend_handle),'uq_Kriging_eval_F')
    % The built-in function handle for evaluating the information matrix F
    % will be used
    
    [trend, Options] = uq_process_option(Options, 'Trend',TrendDefaults, 'struct');
    if trend.Invalid
    error('Invalid trend definition!')
    end
    Trend = trend.Value ;
    
    if trend.Missing
        % The user did not set any options regarding the Trend
        EVT.Type = 'D';
        EVT.Message = 'The default trend options are used.';
        EVT.eventID = 'uqlab:metamodel:kriging:init:trend_defaultsub';
        uq_logEvent(current_model, EVT);
    else
        % Some options regarding the Trend have been set by the user
        % The trendtype MUST be set in order to further define other trend options
        if ~isfield(Trend,'Type') || isempty(Trend.Type)
            EVT.Type = 'W';
            EVT.Message = ['No trend type was defined, therefore the rest of ', ...
                'the options inside the Trend struct are ignored. ', ... 
                'The default trend values are used instead.'];
            EVT.eventID = 'uqlab:metamodel:kriging:init:trend_override';
            uq_logEvent(current_model, EVT);

            Trend = TrendDefaults ;
        end
        if ~ischar(Trend.Type)
            error('Invalid trend type!')
        end

    end

    if isfield(Options, 'Trend')
        % Check for leftover options inside Options.Trend
        uq_options_remainder(Options.Trend, ' Kriging trend options.');
        % Remove Options.Trend
        Options = rmfield(Options,'Trend');
    end

    if INPUT_EXISTS
        [ Trend, EVTs ] = uq_Kriging_initializeTrend( Trend, M, current_model.Internal.Input) ;
    else
        [ Trend, EVTs ] = uq_Kriging_initializeTrend( Trend, M) ;
    end
    % log the returned events, if any
    for iev = 1 : length(EVTs)
        uq_logEvent(current_model, EVTs(iev));
    end
    % Store the trend settings
    current_model.Internal.Kriging.Trend = Trend;
else
    % A user-specified handle for evaluating the information matrix F
    % will be used so by-pass all the additional checks and accept all the
    % supplied options "as is"
    EVT.Type = 'N';
    EVT.Message = sprintf('Trend: using the user-defined function handle: %s',...
        char(trend_handle));
    EVT.eventID = 'uqlab:metamodel:kriging:init:Fhandle_custom';
    uq_logEvent(current_model, EVT);
    
    current_model.Internal.Kriging.Trend = Options.Trend;
    Options = rmfield(Options,'Trend');
    % For now a trend type field is assigned so that this case is easily
    % compatible with the built-in one
    current_model.Internal.Kriging.Trend.Type = 'unknown';
    %% IMPORTANT: If the scaling was not set by the user revert the default to zero
    % to avoid confusion
    if scale.Missing
        EVT.Type = 'W';
        EVT.Message = 'The Scaling option was reverted to 0 because a custom trend is used.';
        EVT.eventID = 'uqlab:metamodel:kriging:init:scaling_revert_custom_trend';
        uq_logEvent(current_model, EVT);
    end
end



%% Parse Correlation function-Related Options
if isfield(Options, 'Corr')
    % Some Correlation related options have been set by the user so
    % parse them
    
    % check whether a non-default function handle is selected for
    % evaluating R:
    [evalRhandle, Options.Corr] = ...
        uq_process_option(Options.Corr, 'Handle',...
        CorrDefaults.Handle, 'function_handle');
    
    % just silently assign the default handle if the user has not specified
    % something
    if evalRhandle.Missing
        % do nothing
    end
    
    if evalRhandle.Invalid
        error('Invalid definition of Correlation function handle!')
    end
    % The rest of the options are only relevant to the default
    % evalR-handle. So only parse them in case the default handle is used:
    if strcmp(char(evalRhandle.Value),'uq_eval_Kernel')
        % first set the handle option
        current_model.Internal.Kriging.GP.Corr.Handle = evalRhandle.Value ;
        
        % Correlation function *type*
        [rtype, Options.Corr] = uq_process_option(Options.Corr, 'Type',...
            CorrDefaults.Type, 'char');
        if rtype.Invalid
            error('Invalid definition of correlation function type!')
        end
        
        current_model.Internal.Kriging.GP.Corr.Type = rtype.Value ;
        
        if rtype.Missing
            msg = sprintf('Correlation function type was set to : %s', ....
                current_model.Internal.Kriging.GP.Corr.Type) ;
            EVT.Type = 'D';
            EVT.Message = msg;
            EVT.eventID = 'uqlab:metamodel:kriging:init:corrfuntype_defaultsub';
            uq_logEvent(current_model, EVT);
            
        end
        
        % Correlation function *family* (it can be either a string for using
        % the built-in ones or a function handle for using a user-defined
        % one
        [rfamily, Options.Corr] = uq_process_option(Options.Corr, 'Family',...
            CorrDefaults.Family, {'char','function_handle'});
        if rfamily.Invalid
            error('Invalid definition of correlation function family!')
        end
        
        current_model.Internal.Kriging.GP.Corr.Family = rfamily.Value ;
        
        if rfamily.Missing
            if strcmpi(class(rfamily.Value),'function_handle')
                msg = sprintf('Correlation family was set to : %s', ....
                    func2str(current_model.Internal.Kriging.GP.Corr.Family)) ;
            else
                msg = sprintf('Correlation family was set to : %s', ....
                    current_model.Internal.Kriging.GP.Corr.Family) ;
            end
            EVT.Type = 'D';
            EVT.Message = msg;
            EVT.eventID = 'uqlab:metamodel:kriging:init:corrfamtype_defaultsub';
            uq_logEvent(current_model, EVT);
        end
        
        % Isotropic
        [risotropic, Options.Corr] = uq_process_option(Options.Corr, 'Isotropic',...
            CorrDefaults.Isotropic, {'double','logical'});
        if risotropic.Invalid
            error('Invalid definition of correlation function''s Isotropic option!')
        end
        
        current_model.Internal.Kriging.GP.Corr.Isotropic = logical(risotropic.Value) ;
        
        if risotropic.Missing
            if current_model.Internal.Kriging.GP.Corr.Isotropic
                msg = sprintf('Correlation function is set to *Isotropic* (default)');
            else
                msg = sprintf('Correlation function is set to *Anisotropic* (default)');
            end
            EVT.Type = 'D';
            EVT.Message = msg;
            EVT.eventID = 'uqlab:metamodel:kriging:init:corrisotropy_defaultsub';
            uq_logEvent(current_model, EVT);
        end
        
        % Nugget
        [nuggetopts, Options.Corr] = uq_process_option(Options.Corr, 'Nugget',...
            CorrDefaults.Nugget, {'double','struct'});
        if nuggetopts.Invalid
            error('Invalid Nugget definition!')
        end
        
        current_model.Internal.Kriging.GP.Corr.Nugget = nuggetopts.Value;
        % No warning message is printed if the user has not set some nugget value
        
        
        % Check for leftover options inside Options.Corr
        uq_options_remainder(Options.Corr, ...
            ' Kriging Correlation function options(.Corr field).');
        
    else
        % If some non-default evalR handle is used, treat all options that
        % are set within the Corr structure as correct and store them
        %
        EVT.Type = 'N';
        EVT.Message = sprintf('Using the user-defined function handle: %s',...
            char(evalRhandle.Value));
        EVT.eventID = 'uqlab:metamodel:kriging:init:corrhandle_custom';
        uq_logEvent(current_model, EVT);
        
        % all the options that were set by the user inside .Corr are stored
        current_model.Internal.Kriging.GP.Corr = Options.Corr; 
        % make sure that the handle option is there
        current_model.Internal.Kriging.GP.Corr.Handle = evalRhandle.Value ;
    end
    % Remove Options.Corr
    Options = rmfield(Options,'Corr');
    
else
    % Default substitution of all options related to the correlation
    % function
    msg = sprintf('The default correlation function options are used:\n%s',...
        printfields(CorrDefaults));
    EVT.Type = 'D';
    EVT.Message = msg;
    EVT.eventID = 'uqlab:metamodel:kriging:init:corrfun_defaultsub';
    uq_logEvent(current_model, EVT);
    
    % No Corr options have been selected so set the default values
    current_model.Internal.Kriging.GP.Corr = CorrDefaults ;
    
end

%% Parse  Estimation Method - related options
if isfield(Options, 'EstimMethod')
    % The Estimation method has been set by the user so parse the related
    % options
    [estmethod, Options] = uq_process_option(Options, 'EstimMethod',...
        EstMethodDefaults.EstimMethod, 'char');
    if estmethod.Invalid
        error('Invalid Hyperparameter estimation method!')
    end
    
    current_model.Internal.Kriging.GP.EstimMethod = estmethod.Value ;
    
    if estmethod.Missing
        msg = sprintf('Hyperparameters estimation method was set to : %s', ...
            current_model.Internal.Kriging.GP.EstimMethod) ;
        EVT.Type = 'D';
        EVT.Message = msg;
        EVT.eventID = 'uqlab:metamodel:kriging:init:estmethod_defaultsub';
        uq_logEvent(current_model, EVT);
    end
  
    switch lower(current_model.Internal.Kriging.GP.EstimMethod)
        % If method is CV make sure that K (is properly defined)
        case 'cv'
            [cvopts, Options] = uq_process_option(Options, 'CV',...
                EstMethodDefaults.CV, 'struct');
            
            if cvopts.Invalid
                error('Invalid Cross-Validation method options!')
            end
            current_model.Internal.Kriging.GP.CV = cvopts.Value ;
            % Make sure that LeaveKOut <= NSamples
            current_model.Internal.Kriging.GP.CV.LeaveKOut = min(current_model.Internal.Kriging.GP.CV.LeaveKOut,...
                current_model.ExpDesign.NSamples);
            
            if cvopts.Missing
                msg = sprintf('Using Cross-Validation method, with: Leave-%i-Out (default).',...
                    current_model.Internal.Kriging.GP.CV.LeaveKOut);
                EVT.Type = 'D';
                EVT.Message = msg;
                EVT.eventID = 'uqlab:metamodel:kriging:init:estmethod_lko_defaultsub';
                uq_logEvent(current_model, EVT);
            end
            
            
    end
else
    % No Estimation method has been set by the user so use the default
    % values
    current_model.Internal.Kriging.GP = merge_structures( ...
        current_model.Internal.Kriging.GP, EstMethodDefaults);
end
% For CV Estimation method calculate the number of classes that corresponds to the
% Leave-K-Out value
if isfield(current_model.Internal.Kriging.GP, 'CV')
    current_model.Internal.Kriging.GP.CV.CV_K = floor(current_model.ExpDesign.NSamples/...
        current_model.Internal.Kriging.GP.CV.LeaveKOut);
end

%% Parse Optimization-Related Options
if isfield(Options, 'Optim')
    % Optimization method
    [optMethod, Options.Optim] = uq_process_option(Options.Optim, ...
        'Method',OptimDefaults.Method, 'char');
    if optMethod.Invalid
        msg = sprintf('Invalid Optimization method. Using the default: %s instead.',...
            optMethod.Value);
        EVT.Type = 'W';
        EVT.Message = msg;
        EVT.eventID = 'uqlab:metamodel:kriging:init:optmethod_override';
        uq_logEvent(current_model, EVT);
        
        current_model.Internal.Kriging.Optim.Method = optMethod.Value ;
    elseif optMethod.Missing
        % Check that the available toolboxes (optim toolbox and global
        % optim toolbox) exist:
        %If global does not exist replace HGA by HCMAES
        % If global and local are not available, replace HGA by CMAES
        try
            % Make sure that the optimization toolbox is available
            evalc('x = fmincon(@(x)x.^2, 0.5, 1, 3);');
            optimization_check = true;
        catch
            optimization_check = false;
        end
        try
            % Make sure that the global optimization toolbox is available
            GAoptions = gaoptimset;
            goptimization_check = true;
        catch
            goptimization_check = false;
        end
        % Should I add a warning in case the default is modified
        if ~optimization_check  && goptimization_check
            optMethod.Value = 'GA' ;
        elseif optimization_check  && ~goptimization_check
            optMethod.Value = 'HCMAES' ;
        elseif ~optimization_check  && ~goptimization_check
            optMethod.Value = 'CMAES' ;
        end
        msg = sprintf('Using the default optimization method: %s.',...
            optMethod.Value);
        EVT.Type = 'D';
        EVT.Message = msg;
        EVT.eventID = 'uqlab:metamodel:kriging:init:optmethod_defaultsub';
        uq_logEvent(current_model, EVT);
        current_model.Internal.Kriging.Optim.Method = optMethod.Value ;
    else
        %make sure that the selected method exists
        if any(strcmpi(OptimKnownMethods, optMethod.Value))
            
            % If the selected method is known use it
            if strcmpi(optMethod.Value, 'lbfgs')
               optMethod.Value = 'BFGS'; 
            end
            % Make sure that corresponding toolbox license is available for
            % the selected optimization method
            try
                % Make sure that the optimization toolbox is avaialble
                evalc('x = fmincon(@(x)x.^2, 0.5, 1, 3);');
                optimization_check = true;
            catch
                optimization_check = false;
            end
            try
                % Make sure that the global optimization toolbox is avaialble
                GAoptions = gaoptimset;
                goptimization_check = true;
            catch
                goptimization_check = false;
            end
            if any(strcmpi(OptimKnownOptimToolboxMethods, optMethod.Value))
                % 'BFGS', 'LBFGS' and 'HSADE' all rely on the 'fmincon function which belongs to the Optimization toolbox
                if ~optimization_check
                    fprintf('The algorithm selected to calibrate the Kriging model is not available\n');
                    fprintf('%s requires the Optimization toolbox which is not available\n',optMethod.Value) ;
                    fprintf('Please select another algorithm or run custom Kriging\n') ;
                    error('Kriging initialization failed: No license for Optimization toolbox') ;
                end
            elseif any(strcmpi(OptimKnownGlobalOptimToolboxMethods,optMethod.Value))
                % 'GA' relies on the ga function which belongs to the global optimization toolbox
                if ~goptimization_check
                    fprintf('The algorithm selected to calibrate the Kriging model is not available\n');
                    fprintf('%s requires the Global Optimization toolbox which is not available\n',optMethod.Value) ;
                    fprintf('Please select another algorithm or run custom Kriging\n') ;
                    error('Kriging initialization failed: No license for Global Optimization toolbox') ;
                end
            elseif any(strcmpi(OptimKnownOptimAndGlobalToolboxMethods,optMethod.Value))
                % 'HGA' relies on ga and fmincon which belong respectively
                % to the Global Optimization toolbox and to the
                % Optimization toolbox
                if ~(optimization_check && goptimization_check)
                    toolbox_result = {'Available', '*Not Available*'};
                    fprintf('The algorithm selected to calibrate the Kriging model is not available\n');
                    fprintf('\t Optimization toolbox: \t\t\t[%s]\n', toolbox_result{2-optimization_check}) ;
                    fprintf('\t Global Optimization toolbox: \t\t[%s]\n', toolbox_result{2-goptimization_check}) ;
                    fprintf('Please select another algorithm or run custom Kriging\n') ;
                    error('Kriging initialization failed: No license for either Optimization toolbox or Global Optimization toolbox') ;             
                end
            else
                % do nothing. Should fall here only if 'sade' is chosen or if the
                % 'OptimKnownMethods is extended
            end
                current_model.Internal.Kriging.Optim.Method = optMethod.Value ;
            
        else
            % If the selected method is unknown raise a warning and use the default
            msg = sprintf('Invalid Optimization method. Using the default: %s instead.',...
                OptimDefaults.Method);
            EVT.Type = 'W';
            EVT.Message = msg;
            EVT.eventID = 'uqlab:metamodel:kriging:init:optmethod_override';
            uq_logEvent(current_model, EVT);

            current_model.Internal.Kriging.Optim.Method = ...
                OptimDefaults.Method ;
        end
    end
    % Get the required fields for the selected Optimization method
    optreqFields = OptimReqFields.(...
        upper(current_model.Internal.Kriging.Optim.Method)) ;
    optmethreqFields = OptimDefaults_Method.(...
        upper(current_model.Internal.Kriging.Optim.Method));
    
    if ~isempty(optmethreqFields)
        optmethreqFields = fieldnames(optmethreqFields) ;
    else
        optmethreqFields = [];
    end
    
    % Try to parse each of the required fields
    for ii = 1 : length(optreqFields)
        [fieldval, Options.Optim] = uq_process_option(Options.Optim, ...
            optreqFields{ii},OptimDefaults.(optreqFields{ii}), ...
            OptimFieldsDatatypes.(upper(optreqFields{ii})));
        switch lower(class(fieldval.Value))
            case 'char'
                printval = fieldval.Value ;
            case 'double'
                if iscolumn(fieldval.Value)
                    printval = uq_sprintf_mat(fieldval.Value') ;
                elseif isrow(fieldval.Value)
                    printval = uq_sprintf_mat(fieldval.Value) ;
                else
                    printval = ['\n', uq_sprintf_mat(fieldval.Value)] ;
                end
            case 'logical'
                printval = num2str(fieldval.Value) ;
            otherwise
                printval = '<not printed>';
        end
        if fieldval.Invalid
            msg = sprintf('Invalid value set at: .Optim.%s! Using the default instead: %s\n',...
                optreqFields{ii}, printval);
            EVT.Type = 'W';
            EVT.Message = msg;
            EVT.eventID = sprintf('uqlab:metamodel:kriging:init:optoptions_override_%i', ...
                ii);
            uq_logEvent(current_model, EVT);
        end
        if fieldval.Missing
            msg = sprintf('The default value for .Optim.%s is used: %s\n',...
                optreqFields{ii}, printval);
            EVT.Type = 'D';
            EVT.Message = msg;
            EVT.eventID = sprintf('uqlab:metamodel:kriging:init:optoptions_defaultsub_%i', ...
                ii);
            uq_logEvent(current_model, EVT);
        end
       
        % set the value
        current_model.Internal.Kriging.Optim.(optreqFields{ii}) = ...
            fieldval.Value ;
    end
    
    % Try to parse optimization method-specific options
    optMethod = upper(current_model.Internal.Kriging.Optim.Method);
    if isfield(Options.Optim, upper(optMethod))
        
        % Some options have been set by the user
        for ii = 1 : length(optmethreqFields)
            [fieldval, Options.Optim.(upper(optMethod))] = uq_process_option(Options.Optim.(upper(optMethod)), ...
                optmethreqFields{ii},...
                OptimDefaults_Method.(upper(optMethod)).(optmethreqFields{ii}), ...
                OptimFieldsDatatypes.(upper(optmethreqFields{ii})));
            switch lower(class(fieldval.Value))
                case 'char'
                    printval = fieldval.Value ;
                case 'double'
                    printval = uq_sprintf_mat(fieldval.Value, '%i') ;
                case 'logical'
                    if fieldval.Value
                        printval = 'true';
                    else
                        printval = 'false';
                    end
                case 'function_handle'
                    printval = func2str(fieldval.Value);
                otherwise
                    printval = '<not printed>';
            end
            if fieldval.Invalid
                msg = sprintf('Invalid value set at: .Optim.%s.%s! Using the default instead: %s\n',...
                    optMethod, optmethreqFields{ii}, printval);
                EVT.Type = 'W';
                EVT.Message = msg;
                EVT.eventID = sprintf('uqlab:metamodel:kriging:init:optmethoptions_override_%i', ...
                    ii);
                uq_logEvent(current_model, EVT);
            end
            if fieldval.Missing
                msg = sprintf('The default value for .Optim.%s.%s is used: %s\n',...
                    optMethod, optmethreqFields{ii}, printval);
                EVT.Type = 'D';
                EVT.Message = msg;
                EVT.eventID = sprintf('uqlab:metamodel:kriging:init:optmethoptions_defaultsub_%i', ...
                    ii);
                uq_logEvent(current_model, EVT);
            end
            % set the value
            current_model.Internal.Kriging.Optim.(upper(optMethod)). ...
                (optmethreqFields{ii}) = fieldval.Value ;
        end
        
        % Check for leftover options inside Options.Optim.(upper(optMethod))
        uq_options_remainder(Options.Optim.(upper(optMethod)), ...
            sprintf(' Kriging Optim.%s options.', optMethod));
        % Remove Options.Optim.(upper(optMethod))
        Options.Optim = rmfield(Options.Optim, optMethod);
        
    else % No Options.Optim have been set by the user
        % For optMethod = 'none' nothing should happen here
        if ~isempty(OptimDefaults_Method.(...
                upper(current_model.Internal.Kriging.Optim.Method)))
            msg = sprintf('The default values for .Optim.%s are used:\n',...
                upper(current_model.Internal.Kriging.Optim.Method));
            
            methdefaults = fieldnames(OptimDefaults_Method.(...
                upper(current_model.Internal.Kriging.Optim.Method)));
            
            for jj = 1 : length(methdefaults)
                switch class(OptimDefaults_Method.(upper(optMethod)). ...
                        (methdefaults{jj}))
                    case 'double'
                        msg = [msg, ...
                            sprintf(' %s : %s\n', methdefaults{jj},...
                            num2str(OptimDefaults_Method.(upper(optMethod)). ...
                            (methdefaults{jj})))];
                    case 'char'
                        msg = [msg, ...
                        sprintf(' %s : %s\n', methdefaults{jj},...
                            OptimDefaults_Method.(upper(optMethod)). ...
                            (methdefaults{jj}))];
                    otherwise
                        msg = [msg, ...
                        sprintf(' %s : %s\n', methdefaults{jj},...
                            ['<',class(methdefaults{jj}),'>'])];
                end
            end
            % Log the default substitution event
            EVT.Type = 'D';
            EVT.Message = msg;
            EVT.eventID = sprintf('uqlab:metamodel:kriging:init:optimsome_defaultsub');
            uq_logEvent(current_model, EVT);
            
            % No options have been set. Set the defaults
            current_model.Internal.Kriging.Optim.(...
                upper(current_model.Internal.Kriging.Optim.Method))...
                = OptimDefaults_Method.(...
                upper(current_model.Internal.Kriging.Optim.Method));
        end
    end
    
    
    % Additional checks for some specific fields
    current_model = validateOptimOptions(current_model);
    if isfield(Options.Optim, 'InitialValue')
        Options.Optim = rmfield(Options.Optim, 'InitialValue');
    end
    if isfield(Options.Optim, 'Bounds')
        Options.Optim = rmfield(Options.Optim, 'Bounds');
    end
    
    
    % Check for leftover options inside Options.Optim
    uq_options_remainder(Options.Optim, ' Kriging Optim options.');
    % Remove Options.Optim
    Options = rmfield(Options, 'Optim');
    
else
    
        % Check that the available toolboxes (optim toolbox and global
        % optim toolbox) exist:
        % If global does not exist replace HGA by HCMAES
        % If global and local are not available, replace HGA by CMAES
        try
            % Make sure that the optimization toolbox is avaialble
            evalc('x = fmincon(@(x)x.^2, 0.5, 1, 3);');
            optimization_check = true;
        catch
            optimization_check = false;
        end
        try
            % Make sure that the global optimization toolbox is avaialble
            GAoptions = gaoptimset;
            goptimization_check = true;
        catch
            goptimization_check = false;
        end
        % Should I add a warning in case the default is modified
        if ~optimization_check  && goptimization_check
            OptimDefaults.Method = 'GA' ;
        elseif optimization_check  && ~goptimization_check
            OptimDefaults.Method = 'HCMAES' ;
        elseif ~optimization_check  && ~goptimization_check
            OptimDefaults.Method = 'CMAES' ;
        end
        
    current_model.Internal.Kriging.Optim = OptimDefaults ;
    % Additional checks for some specific fields
    current_model = validateOptimOptions(current_model);
    % No Optimization options have been selected so set the defaults
    msg = sprintf('The default values for .Optim are used:\n%s', ...
        printfields(current_model.Internal.Kriging.Optim));
    EVT.Type = 'D';
    EVT.Message = msg;
    EVT.eventID = sprintf('uqlab:metamodel:kriging:init:optimall_defaultsub');
    uq_logEvent(current_model, EVT);
end

%% Check for unused options
% Remove some fields that are not processed here:
fieldsToRemove = skipFields(isfield(Options,skipFields)) ;

Options = rmfield(Options, fieldsToRemove);

% Check if there was something else provided:
uq_options_remainder(Options, ...
    current_model.Name, ...
    skipFields, current_model);


%% Add the property where the main Kriging results are stored
uq_addprop(current_model, 'Kriging');
success = 1;

end % END OF uq_Kriging_initialize FUNCTION



%% Helper functions for the Kriging initialization
% (To be moved elsewhere)

function sout = merge_structures(varargin)
%MERGE_STRUCTURES merges structure variables given that there is no
%overlap on their fields

% collect the names of the fields
fnames = [];
for k = 1:nargin
    try
        fnames = [fnames ; fieldnames(varargin{k})];
    catch 
        % do nothing
    end
end

% Make sure the field names are unique
if numel(fnames) ~= numel(unique(fnames))
    error('Internal Kriging initialization error: Field names must be unique!');
end

% Now concatenate the data from each structure into a cell array
cellarr = {};
for k = 1:nargin
    cellarr = [cellarr ; struct2cell(varargin{k})];
end

% transform the concatenated data from cell to struct
sout = cell2struct(cellarr, fnames, 1);
end


function msg = printfields(S, tabs)
% PRINTFIELDS prints the names of the fields of a structure depending on their data
% type
if nargin < 2
    tabs = '';
end
msg = '';
fnames = fieldnames(S);
for jj = 1 : length(fnames)
    switch class(S.(fnames{jj}))
        case 'char'
            msg_new = sprintf('%s %s : %s\n', tabs, fnames{jj}, ...
                S.(fnames{jj})) ;
        case 'double'
            msg_new = sprintf('%s %s : %s\n', tabs, fnames{jj}, ...
                uq_sprintf_mat(S.(fnames{jj}))) ;
        case 'struct'
            msg_new = sprintf('%s %s(contents) : \n', tabs, fnames{jj});
            tabs = sprintf('%s\t',tabs);
            msg_new = [msg_new, printfields(S.(fnames{jj}), tabs)] ;
            tabs = '';
        otherwise
            msg_new = sprintf('%s %s : %s\n', tabs, fnames{jj}, ...
                '<not printed>') ;
    end
    msg = [msg, msg_new];
end

end

function current_model = validateOptimOptions(current_model)
% helper function to validate and fix the dimension of the bounds or
% initial value
% It currently handles 3 possible cases:
% 1) Using built-in uq_eval_Kernel and anisotropic -> theta is Mx1
% vector,  bounds is 2xM 
% 2) Using built-in uq_eval_Kernel and isotropic -> theta is 1x1
% scalar,  bounds is 2x1 
% 3) Using user-defined eval_R function -> no checks on theta and bounds!

if strcmp(char(current_model.Internal.Kriging.GP.Corr.Handle),...
        'uq_eval_Kernel')
    M = current_model.Internal.Runtime.M;
    %% Initial Value
    if isfield(current_model.Internal.Kriging.Optim, 'InitialValue' )
        th0 = current_model.Internal.Kriging.Optim.InitialValue ;
        % get the Isotropic flag of the correlation function
        isIsotropic = current_model.Internal.Kriging.GP.Corr.Isotropic;
        
        % Make sure that Optim.InitialValue has proper dimensions
        if isIsotropic
            dimensionCheckOK = isscalar(th0);
            assert(dimensionCheckOK,...
                'For isotropic correlation function the hyperparameter defined in .Optim.InitialValue is expected to be a scalar!');
        else
            dimensionCheckOK = sum(size(th0) == [M,1])== 2 ||...
                sum( size(th0) == [1,M])== 2;
            if ~dimensionCheckOK
                % if a scalar theta0 is assigned and M>1, replicate theta0 across
                % all dimensions
                if isscalar(th0)
                    current_model.Internal.Kriging.Optim.InitialValue = repmat(th0,M,1) ;
                    % Print the InitialValue that is going to be used
                    msg = sprintf('\t> Optimization initial value is updated to:\n%s\n',...
                        uq_sprintf_mat(current_model.Internal.Kriging.Optim.InitialValue));
                    EVT.Type = 'N';
                    EVT.Message = msg;
                    EVT.eventID = 'uqlab:metamodel:kriging:init:optinitval_replicate';
                    uq_logEvent(current_model, EVT);
                else
                    % theta0 is neither scalar or M dimensional vector, so raise an
                    % error
                    error('Dimension mismatch between Optim.InitialValue and Experimental Design!')
                end
            end
        end
        
        
    end
    
    %% Bounds
    if isfield(current_model.Internal.Kriging.Optim, 'Bounds' )
        % Make sure that Optim.Bounds have proper dimensions
        thLB  = current_model.Internal.Kriging.Optim.Bounds(1,:);
        thUB  = current_model.Internal.Kriging.Optim.Bounds(2,:);
        % get the Isotropic flag of the correlation function
        isIsotropic = current_model.Internal.Kriging.GP.Corr.Isotropic;
        % a flag that determines whether some update on the bounds occured during validation. 
        % If this true then the updated bounds are printed on the screen         
        printoutUpdate = false;
        
        if isIsotropic
            dimensionCheckOK = isscalar(thUB) & isscalar(thLB);
            assert(dimensionCheckOK,...
                'For isotropic correlation function .Optim.Bounds are expected to be a 2x1 vector!');
            
        else
            dimensionCheckUBOK = sum(size(thUB) == [M,1])==2 ||...
                sum( size(thUB) == [1,M])==2;
            dimensionCheckLBOK = sum(size(thLB) == [M,1])==2 ||...
                sum( size(thLB) == [1,M])==2;
            
            if ~dimensionCheckUBOK
                % if a scalar thUB is assigned and M>1, replicate thUB across
                % all dimensions
                if sum(size(thUB) == [1,1]) == 2
                    thUB = repmat(thUB,1,M) ;
                    printoutUpdate = true;
                else
                    % thUB is neither scalar or M dimensional vector, so raise an
                    % error
                    error('Dimension mismatch between Optim.Bounds and Experimental Design!')
                end
            end
            
            if ~dimensionCheckLBOK
                % if a scalar thLB is assigned and M>1, replicate thLB across
                % all dimensions
                if sum(size(thLB) == [1,1]) == 2
                    thLB = repmat(thLB,1,M) ;
                    printoutUpdate = true;
                else
                    % thLB is neither scalar or M dimensional vector, so raise an
                    % error
                    error('Dimension mismatch between Optim.Bounds and Experimental Design!')
                end
            end
        end
        current_model.Internal.Kriging.Optim.Bounds = [thLB ; thUB] ;
        if printoutUpdate
            % Print the Bounds that are going to be used
            msg = sprintf('\t> Optimization bounds value is updated to:\n%s\n',...
                uq_sprintf_mat(current_model.Internal.Kriging.Optim.Bounds));
            EVT.Type = 'N';
            EVT.Message = msg;
            EVT.eventID = 'uqlab:metamodel:kriging:init:optbounds_replicate';
            uq_logEvent(current_model, EVT);
        end
    end
end

end




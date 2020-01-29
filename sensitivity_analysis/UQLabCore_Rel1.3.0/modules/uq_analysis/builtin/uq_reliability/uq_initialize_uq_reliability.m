function success = uq_initialize_uq_reliability(current_analysis)
% success = UQ_INITIALIZE_UQ_RELIABILITY(current_analysis):
%     initializes the structural reliability analyses defined in 
%     CURRENT_ANALYSIS. All initializations are done in this file.
% 
% See also: UQ_RELIABILITY

success = 0;


%%
% User options
Options = current_analysis.Options; 

% Actual options that the module will use
Internal = current_analysis.Internal; 

%% Check if the user has provided a Model and / or an Input
[opt, Options] = uq_process_option(Options, 'Model', uq_getModel, 'uq_model');
Internal.Model = opt.Value;

[opt, Options] = uq_process_option(Options, 'Input', uq_getInput, 'uq_input');
Internal.Input = opt.Value;

%% Set the Method:
% Depending on the method, some parts of the initialization script are
% needed or not:
InitSimulation = false;
InitFORM = false;
InitGradient = false;
InitIS = false;
InitSubsetSim = false;
InitAKMCS = false;
InitInvFORM = false ;
%% Reset the limitState count
uq_evalLimitState('reset');

%% Set the method specify defaults
%process the method
[Method, Options] = uq_process_option(Options, 'Method', 'mc', 'char');

if Method.Missing || Method.Invalid
    error('\nAn analysis method should be provided.\n');
else
    % Map the multiple methods description to a fixed set:
    % mc, is, form, sorm
    switch lower(Method.Value)
        case {'mc','mcs','montecarlo', 'monte carlo'}
            Internal.Method = 'mc';
            % Needs to initialize the simulation options:
            InitSimulation = true;
            DEFAULTSaveEvaluations = 1;
            DEFAULTBatchSize = 1e4;
            
        case {'importancesampling','is','importance sampling'}
            Internal.Method = 'is';
            % Needs to initialize the simulation options, FORM and the
            % gradient:
            InitSimulation = true;
            InitFORM = true;
            InitGradient = true;
            InitIS = true;
            DEFAULTSaveEvaluations = 1;
            DEFAULTBatchSize = 1e2;
            
        case 'form'
            [opt, Options] = uq_process_option(Options, 'NoSORM');
            if opt.Disabled || opt.Missing || opt.Invalid
                Internal.NoSORM = true;
            else
                Internal.NoSORM = opt.Value;
            end
            
            % Decide whether to change to SORM or not:
            if Internal.NoSORM
                Internal.Method = 'form';
            else
                Internal.Method = 'sorm';
            end
            % Needs to initialize FORM and the gradient:
            InitFORM = true;
            InitGradient = true;
            DEFAULTSaveEvaluations = 1;
            
        case 'sorm'
            Internal.Method = 'sorm';
            % Needs to initialize FORM and the gradient:
            InitFORM = true;
            InitGradient = true;
            DEFAULTSaveEvaluations = 1;
            
        case 'subset'
            % Needs to initialize simulation options and subset simulation
            % options
            Internal.Method = 'subset';
            InitSimulation = true;
            InitSubsetSim = true;
            DEFAULTSaveEvaluations = 1;
            DEFAULTBatchSize = 1e3;
            
        case 'akmcs'
            % Needs to initialize simulation options and ak-mcs options
            Internal.Method = 'akmcs';
            InitSimulation = true;
            InitAKMCS = true;
            DEFAULTSaveEvaluations = 1;
            DEFAULTBatchSize = 1e4;
            DEFAULTNugget = 1e-10; 
        case 'apckmcs'
            % Needs to initialize simulation options and ak-mcs options
            Internal.Method = 'akmcs';
            Options.AKMCS.MetaModel = 'PCK';
            InitSimulation = true;
            InitAKMCS = true;
            DEFAULTSaveEvaluations = 1;
            DEFAULTBatchSize = 1e4;
            DEFAULTNugget = 1e-10; 
            
        case { 'inverseform','iform' }
            [opt, Options] = uq_process_option(Options, 'NoSORM');
            if opt.Disabled || opt.Missing || opt.Invalid
                Internal.NoSORM = true;
            else
                Internal.NoSORM = opt.Value;
            end
            
            % Decide whether to change to SORM or not:
            if Internal.NoSORM
                Internal.Method = 'inverseform';
            else
                Internal.Method = 'sorm';
            end
            % Needs to initialize FORM and the gradient:
            InitFORM = true;
            InitGradient = true;
            InitInvFORM = true ;
            DEFAULTSaveEvaluations = 1;
        otherwise
            error('\nThe provided analysis method "%s" is not correct.', ...
                Method.Value);
            
    end
end





%% COMMON OPTIONS:

% Process the limit state function (as struct):
DefLS.Threshold = 0;
DefLS.CompOp = '<=';

[opt, Options] = uq_process_option(Options, 'LimitState', DefLS, 'struct');
LimitState = opt.Value;

% Process the fields of LimitState:
% - Threshold:
[opt, LimitState] = uq_process_option(LimitState, 'Threshold', 0, 'double');
LimitState.Threshold = opt.Value;

% - Comparison option (CompOp):
[opt, LimitState] = uq_process_option(LimitState, 'CompOp', '<=', 'char');
switch opt.Value
    case {'<', '<=', 'leq', '>', '>=', 'geq'}
        LimitState.CompOp = opt.Value;
    otherwise
        LimitState.CompOp = opt.Default;
end


% LimitState finished processing now:
Internal.LimitState = LimitState;

% Save the evaluations (when possible)
[opt, Options] = uq_process_option(Options, 'SaveEvaluations', DEFAULTSaveEvaluations, {'logical', 'double', 'char'});
switch opt.Value
    case {1, true, 'yes','all','enabled', 'debug', 'true'}
        Internal.SaveEvaluations = 1;
    otherwise
        Internal.SaveEvaluations = opt.Default;
end

% Transient analysis (for now, only FORM treats this differently)
[opt, Options] = uq_process_option(Options, 'Transient', false, {'logical', 'double', 'char'});
switch opt.Value
    case {1, true, 'yes','all','enabled'}
        Internal.Transient = 1;
    otherwise
        Internal.Transient = opt.Default;
end

% Initialize the display:
[Options, Internal] = uq_initialize_display(Options, Internal);

%%%UQHPCSTART%%%
[Options, Internal] = uq_initialize_HPC(Options, Internal);
% Set the default granularity options:
% By default, only simulation methods are distributed.
HPCOptions = Internal.HPC;

if HPCOptions.Enabled
    
    [opt, HPCOptions] = uq_process_option(HPCOptions, 'FORM', false, {'double', 'logical'});
    Internal.HPC.FORM = opt.Value;
    
    [opt, HPCOptions] = uq_process_option(HPCOptions, 'SORM', false, {'double', 'logical'});
    Internal.HPC.SORM = opt.Value;
    
    [opt, HPCOptions] = uq_process_option(HPCOptions, 'MC', true, {'double', 'logical'});
    Internal.HPC.MC = opt.Value;
    
    [opt, HPCOptions] = uq_process_option(HPCOptions, 'IS', true, {'double', 'logical'});
    Internal.HPC.IS = opt.Value;
    
    % Remove the option thas was processed before:
    HPCOptions = rmfield(HPCOptions, 'Enabled');
    
    % Check if there was something else provided:
    uq_options_remainder(HPCOptions, ...
        'the HPC configuration of the analysis module');
    
else
    Internal.HPC.FORM = false;
    Internal.HPC.SORM = false;
    Internal.HPC.MC = false;
    Internal.HPC.IS = false;
end
%%%UQHPCEND%%%

OptionsFields = fieldnames(Options);


%% SIMULATION METHODS (MC, IS, subset sim., and AKMCS)

if InitSimulation
    
    SimOpts = struct;
    
    [opt, Options] = uq_process_option(Options, 'Simulation', SimOpts, 'struct');
    SimOpts = opt.Value;
    
    
    % Set the batch size to evaluate at a time
    [opt, SimOpts] = uq_process_option(SimOpts, 'BatchSize', DEFAULTBatchSize, 'double');
    Internal.Simulation.BatchSize = opt.Value;
    
    % Set the sampling method
    [opt, SimOpts] = uq_process_option(SimOpts, 'Sampling', 'mc', 'char');
    Internal.Simulation.Sampling = opt.Value;
    
    % Set the target coefficient of variation
    [opt, SimOpts] = uq_process_option(SimOpts, 'TargetCoV', [], 'double');
    Internal.Simulation.TargetCoV = opt.Value;
    
    % If there is a target CoV, we can set max. runs to Inf by default,
    % otherwise, we limit it to 10^5:
    if isempty(opt.Value)
        if strcmp(current_analysis.Options.Method, 'IS')
            MRdefault = 1e3;
        else
            MRdefault = 1e5;
        end
    else
        MRdefault = inf;
    end
    
    % Maximum number of function evaluations
    [opt, SimOpts] = uq_process_option(SimOpts, 'MaxSampleSize', MRdefault, 'double');
    Internal.Simulation.MaxSampleSize = opt.Value;
    if opt.Invalid || opt.Missing || opt.Disabled
        fprintf('Warning: Maximum number of model evaluations was limited to 10^5\n');
    end
    
    
    % Confidence level Alpha:
    [opt, SimOpts] = uq_process_option(SimOpts, 'Alpha', 0.05, 'double');
    Internal.Simulation.Alpha = opt.Value;
    
    uq_options_remainder(SimOpts, 'the simulation methods');
    
end

%% Options for the FORM algorithm (and therefore also for SORM and IS)
if InitFORM
    
    FORMOpts = struct;
    [opt, Options] = uq_process_option(Options, 'FORM', FORMOpts, 'struct');
    FORMOpts = opt.Value;
    
    
    % FORM algorithm Specific Options
    [opt, FORMOpts] = uq_process_option(FORMOpts, 'StartingPoint', [], 'double');
    Internal.FORM.StartingPoint = opt.Value;

    % Solution algorithm:
    [opt, FORMOpts] = uq_process_option(FORMOpts, 'Algorithm', 'iHLRF', 'char');
    
    % Only two specific strings are accepted:
    if ~any(strcmpi(opt.Value, {'hlrf', 'ihlrf'}))
        Internal.FORM.Algorithm = opt.Default;
        opt.Invalid = true;
        % Show a warning if the user provided an invalid string:
        fprintf('\nWarning: The provided algorithm "%s" was not recognized. Switching to "%s".\n', opt.Value, opt.Default);
    else
        Internal.FORM.Algorithm = opt.Value;
    end
    
    
    % Stop criteria on U: (Stagnation)
    [opt, FORMOpts] = uq_process_option(FORMOpts, 'StopU', 1e-4, 'double');
    Internal.FORM.StopU = opt.Value;
    
    
    % Stop criteria on G: (Convergence)
    [opt, FORMOpts] = uq_process_option(FORMOpts, 'StopG', 1e-6, 'double');
    Internal.FORM.StopG = opt.Value;
    
    
    % Maximum number of iterations:
    [opt, FORMOpts] = uq_process_option(FORMOpts, 'MaxIterations', 1e2, 'double');
    Internal.FORM.MaxIterations = opt.Value;
    
    
end

if InitInvFORM
    invFORMOpts = struct;
    
    [opt, Options] = uq_process_option(Options, 'invFORM', invFORMOpts, 'struct');
    invFORMOpts = opt.Value;

        % inverse FORM algorithm Specific Options
    [opt, invFORMOpts] = uq_process_option(invFORMOpts, 'StartingPoint', [], 'double');
    Internal.invFORM.StartingPoint = opt.Value;
    
    % Solution algorithm:
    [opt, invFORMOpts] = uq_process_option(invFORMOpts, 'Algorithm', 'AMV', 'char');
    % Only two specific strings are accepted:
    if ~any(strcmpi(opt.Value, {'amv', 'cmv', 'hmv','sqp'}))
        Internal.invFORM.Algorithm = opt.Default;
        opt.Invalid = true;
        % Show a warning if the user provided an invalid string:
        fprintf('\nWarning: The provided algorithm "%s" was not recognized. Switching to "%s".\n', opt.Value, opt.Default);
    else
        Internal.invFORM.Algorithm = opt.Value;
    end

    % Target reliability index
    [opt, invFORMOpts] = uq_process_option(invFORMOpts, 'TargetBetaHL', 1.6449, 'double');
    Internal.invFORM.TargetBetaHL = opt.Value;
    
    % Stop criteria on U: (Stagnation)
    [opt, invFORMOpts] = uq_process_option(invFORMOpts, 'StopU', 1e-3, 'double');
    Internal.invFORM.StopU = opt.Value;
    
    % Stop criteria on U: (Stagnation)
    [opt, invFORMOpts] = uq_process_option(invFORMOpts, 'StopG', 1e-3, 'double');
    Internal.invFORM.StopG = opt.Value;
 
    % Maximum number of iterations:
    [opt, invFORMOpts] = uq_process_option(invFORMOpts, 'MaxIterations', 1e2, 'double');
    Internal.invFORM.MaxIterations = opt.Value;
    
end
%% Differentiation options (also for FORM, SORM, IS):
if InitGradient
    % default to standardized calculation of the gradient
    if isfield(Options, 'Gradient')
        GradOpts = Options.Gradient;
        if ~isfield(GradOpts, 'Step') || isempty(GradOpts.Step)
            Options.Gradient.Step = 'Standardized';
        end
    else
        Options.Gradient.Step = 'Standardized';
    end
    [Options, Internal] = uq_initialize_gradient(Options, Internal, true);
end

%% Importance sampling (IS) special options:
if InitIS
    
    % See if the option was provided at all:
    ISOpts.FORM = [];
    ISOpts.Instrumental = [];
    [opt, Options] = uq_process_option(Options, 'IS', ISOpts, 'struct');
    ISOpts = opt.Value;
    
    % Process the field FORM:
    [opt, ISOpts] = uq_process_option(ISOpts, 'FORM', [], {'struct', 'uq_analysis', 'double'});
    Internal.IS.FORM = opt.Value;
    % In case the provided option is a double other than []
    if isnumeric(opt.Value) && ~isempty(opt.Value)
        fprintf('\nWarning: Invalid option for IS.FORM, it will be ignored.\n');
        opt.Value = [];
    end
    
    % If the user gave an analysis, check if it is solved and retrieve the
    % results:
    if isa(opt.Value, 'uq_analysis')
        if isprop(opt.Value, 'Results')
            opt.Value = opt.Value.Results(end);
            Internal.IS.FORM = opt.Value;
        else
            fprintf('\nWarning: The FORM results given to importance sampling are not valid.\n');
            fprintf('FORM will be executed to find the design point.\n');
            Internal.IS.FORM = [];
            opt.Value = [];
        end
    end


    % Check that if it was given by the user, it has the required fields,
    % otherwise run FORM later:
    if isstruct(opt.Value)
        if ~isfield(opt.Value, 'Ustar') % || ~isfield(opt.Value, 'BetaHL')
            fprintf('\nWarning: The FORM results given to importance sampling are not valid.\n');
            fprintf('FORM will be executed to find the design point.\n');
            Internal.IS.FORM = [];
        end
    end
    
    % Process the instrumental density, field SOpts.IS.Instrumental
    [opt, ISOpts] = uq_process_option(ISOpts, 'Instrumental', [], {'struct', 'uq_input', 'double'});
    
    if isnumeric(opt.Value) && ~isempty(opt.Value) % In case the provided option is a double other than []
        fprintf('\nWarning: Invalid option for IS.Instrumental, it will be ignored.\n');
        Internal.IS.Instrumental = [];
        
    elseif length(opt.Value)==1 % There is only one instrumental provided
        Internal.IS.Instrumental = opt.Value;
        if isstruct(opt.Value)
            try
                InstrumentalInput = uq_createInput(opt.Value, '-private');
                Internal.IS.Instrumental = InstrumentalInput;
            catch ME
                fprintf('\nWarning: The struct given in IS.Instrumental does not define a proper instrumental density');
                fprintf('\nInput could not be initialized with error:\n%s\n', ME.message);
                fprintf('\nFORM will be executed to find a proper instrumental density\n');
                Internal.IS.Instrumental = [];
            end
        elseif isa(opt.Value,'uq_input')
            Internal.IS.Instrumental = opt.Value;
        else
            fprintf('\nWarning: Invalid option for IS.Instrumental, it will be ignored.\n');
            Internal.IS.Instrumental = [];
        end
        
    elseif isstruct(opt.Value) || isa(opt.Value,'uq_input') % struct contains multiple InstrInputs or InstrOpts
        
            for uu = 1:length(opt.Value) % turn all Opts into Inputs
                if isa(opt.Value(uu),'uq_input') % already in desired form
                    Internal.IS.Instrumental(uu) = opt.Value(uu);
                elseif isstruct(opt.Value(uu)) % turn Opts into Input
                    try
                        Internal.IS.Instrumental(uu) = uq_createInput(opt.Value(uu), '-private');                        
                    catch ME
                        fprintf('\nWarning: The struct given in IS.Instrumental does not define a proper instrumental density');
                        fprintf('\nInput could not be initialized with error:\n%s\n', ME.message);
                        fprintf('\nFORM will be executed to find a proper instrumental density\n');
                        Internal.IS.Instrumental = [];
                        break
                    end
                else % structure contains nothing useful
                    fprintf('\nWarning: Invalid option for IS.Instrumental, it will be ignored.\n');
                    Internal.IS.Instrumental = [];
                    break
                end
            end
    elseif isempty(opt.Value)
        % No Instrumental distribution provided, do nothing
        Internal.IS.Instrumental = [];
    else
        fprintf('\nWarning: The struct given in IS.Instrumental does not define a proper instrumental density');
        fprintf('\nInput could not be initialized');
        fprintf('\nFORM will be executed to find a proper instrumental density\n');
        Internal.IS.Instrumental = [];
    end
        
    % Check if both were given, and in that case prefer the instrumental
    % density:
    if ~isempty(Internal.IS.Instrumental) && ~isempty(Internal.IS.FORM)
        fprintf('\nWarning: Importance sampling cannot accept both FORM results and an instrumental density\n');
        fprintf('FORM results will be ignored\n');
        Internal.IS.FORM = [];
    end
    % Check if there were more options given:
    uq_process_option(ISOpts, ' importance sampling.');
end

%% Initialize subset simulation options:
% Note that the default values are chosen in correpsondance to: 
%   Au and Beck (2001), Estimation of small failure probabilities in high
%   dimensions by subset simulation, Prob. Eng. Mech. 16(4), p 263-277.
if InitSubsetSim
    
    SSOpts = struct;
    
    [opt, Options] = uq_process_option(Options, 'Subset', SSOpts, 'struct');
    SSOpts = opt.Value;
    
    %set the intermediate failure probability p0
    [opt, SSOpts] = uq_process_option(SSOpts, 'p0', 0.1, 'double');
    Internal.Subset.p0 = opt.Value;
    
    %set the maximum number of subsets
    MaxSubsets = floor(Internal.Simulation.MaxSampleSize / Internal.Simulation.BatchSize / (1-Internal.Subset.p0));
    [opt, SSOpts] = uq_process_option(SSOpts, 'MaxSubsets', MaxSubsets,'double');
    Internal.Subset.MaxSubsets = min(opt.Value, MaxSubsets);
    
    %set the MH acceptance criterion
    [opt, SSOpts] = uq_process_option(SSOpts, 'Componentwise', 1, 'double');
    Internal.Subset.Componentwise = opt.Value;
    
    %set the random walker of the Markov Chain
    RWOpts = struct;
    [opt, SSOpts] = uq_process_option(SSOpts, 'Proposal', RWOpts, 'struct');
    RWOpts = opt.Value;
    
    %set the proposal distribution of the random walker
    [opt, RWOpts] = uq_process_option(RWOpts, 'Type', 'uniform', 'char');
    Internal.Subset.Proposal.Type = opt.Value;
    
    % set the proposal distribution of the random walker
    [opt, RWOpts] = uq_process_option(RWOpts, 'Parameters', 1, 'double');
    Internal.Subset.Proposal.Parameters = opt.Value;
    
end

%% Initialize AK-MCS options
if InitAKMCS
    AKOpts = struct;
    
    [opt, Options] = uq_process_option(Options, 'AKMCS', AKOpts, 'struct');
    AKOpts = opt.Value;
    
    %check the initial experimental design
    if ~isfield(AKOpts, 'IExpDesign')
        AKOpts.IExpDesign.N = max(10, 2*length(Internal.Input.Marginals));
        AKOpts.IExpDesign.Sampling = 'mc';
    end
    if ~isfield(AKOpts.IExpDesign, 'X')
        IEDefault.N = max(10, 2*length(Internal.Input.Marginals));
        IEDefault.Sampling = 'mc';
        [opt, AKOpts] = uq_process_option(AKOpts,'IExpDesign',IEDefault,'struct');
        Internal.AKMCS.IExpDesign = opt.Value;
    else 
        [opt, AKOpts] = uq_process_option(AKOpts,'IExpDesign',[],'struct');
        Internal.AKMCS.IExpDesign = opt.Value;
    end
    
    %check the meta-model type
    [opt, AKOpts] = uq_process_option(AKOpts, 'MetaModel', 'Kriging', 'char');
    MetaModel = opt.Value;
    
    %put the meta-model type name as a standardized on for consistencies
    switch lower(MetaModel)
        case 'kriging'
            MetaModel = 'Kriging';
        case 'pck'
            MetaModel = 'PCK';
    end
    Internal.AKMCS.MetaModel = MetaModel;
    
    %process the Kriging options if available
    [opt, AKOpts] = uq_process_option(AKOpts, MetaModel, [], 'struct');
    Internal.AKMCS.(MetaModel) = opt.Value;
    
    switch lower(MetaModel)
        case 'kriging'
            %check for the Nugget (in Kriging)
            if isfield(Internal.AKMCS.Kriging, 'Corr')
                if ~isfield(Internal.AKMCS.Kriging.Corr, 'Nugget')
                    Internal.AKMCS.Kriging.Corr.Nugget = DEFAULTNugget;
                end
            else
                Internal.AKMCS.Kriging.Corr.Nugget = DEFAULTNugget;
            end
            
        case 'pck'
            %check for the Nugget (in Kriging)
            if isfield(Internal.AKMCS.PCK,'Kriging')&& isfield(Internal.AKMCS.PCK.Kriging, 'Corr')
                if ~isfield(Internal.AKMCS.PCK.Kriging.Corr, 'Nugget')
                    Internal.AKMCS.PCK.Kriging.Corr.Nugget = DEFAULTNugget;
                end
            else
                Internal.AKMCS.PCK.Kriging.Corr.Nugget = DEFAULTNugget;
            end
            
    end
    %check whether there are Kriging options specified
    %-> done directly in uq_akmcs at the moment
    
    %check the learning function
    [opt, AKOpts] = uq_process_option(AKOpts, 'LearningFunction', 'U', 'char');
    Internal.AKMCS.LearningFunction = opt.Value;
    
    %check the number of total added samples
    [opt, AKOpts] = uq_process_option(AKOpts, 'MaxAddedED', 1000, 'double');
    Internal.AKMCS.MaxAddedSamplesTotal = opt.Value;
    
    %locking the number of samples in each batch
    Internal.AKMCS.MaxAddedSamplesInBatch = Internal.AKMCS.MaxAddedSamplesTotal;
    
    %check the convergence criterion
    [opt, AKOpts] = uq_process_option(AKOpts, 'Convergence', 'stopU', 'char');
    Internal.AKMCS.Convergence = opt.Value;
    
    %locking the BatchSize = MaxNumberOfSamples
    Internal.Simulation.BatchSize = Internal.Simulation.MaxSampleSize;

end

%% Remove the type and name and  check the reminder
% Remove the option Type, that is not processed:
Options = rmfield(Options, 'Type');

% Remove also the Name option, if provided:
if isfield(Options, 'Name')
    Options = rmfield(Options, 'Name');
end

% Check if there was something else provided:
uq_options_remainder(Options, ...
    'Reliability');

%% Give the filtered options back to the object
current_analysis.Internal = Internal;
success = 1;

%% Run the analysis
uq_runAnalysis(current_analysis);
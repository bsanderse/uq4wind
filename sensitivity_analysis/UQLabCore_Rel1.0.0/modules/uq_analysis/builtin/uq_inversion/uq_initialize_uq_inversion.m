function [Success] = uq_initialize_uq_inversion(CurrentAnalysis)
% UQ_INITIALIZE_UQ_INVERSION initializes an inverse analysis in UQLab by
%   going through the user specified and default options.
%
%   See also: UQ_INVERSION, UQ_INVERSION_AUGMENTINPUT, 
%   UQ_INVERSION_NORMPDF, UQ_INVERSION_LOGNORMPDF, UQ_INVERSION_LIKELIHOOD 
%   UQ_EVAL_LOGPDF, UQ_INITIALIZE_UQ_INVERSION_PROPOSAL

%% SUCCESS
Success = 0;

%% OPTIONS
% User options
Options = CurrentAnalysis.Options;

% Actual options that the module will use
Internal = CurrentAnalysis.Internal; 


%% CUSTOM LIKELIHOOD SWITCH
[Opt,Options] = uq_process_option(Options,'LogLikelihood');
if Opt.Missing
    % customLikelihood flag false
    Internal.customLikeli = false;
    CUSTOM_LIKELI = false;
else
    % custom log likelihood specified
    Internal.customLikeli = true;
    CUSTOM_LIKELI = true;
        % make sure only one data group is provided
    if length(Options.Data) ~= 1
        error('Multiple data groups not supported for custom logLikelihood')
    end
    
    % Assign number of data groups to internal
    Internal.nDataGroups = length(Options.Data);
    Internal.LogLikelihood = Opt.Value;
end

%% DISPLAY
% Set the verbosity level:
[Options, Internal] = uq_initialize_display(Options, Internal);

%% INPUT
[Opt,Options] = uq_process_option(Options,'Prior');
% missing & invalid
if Opt.Missing
    % check that no input object was specified in the discrepancy options
    if isfield(Options,'Discrepancy')
        if isfield(Options.Discrepancy,'Prior')
            error('The prior needs to be specified, if uq_input object is used in discrepancy options')
        end
    end
    % check that no input object was specified in the sampler options
    if isfield(Options,'Solver')
        if isfield(Options.Solver,'MCMC')
            if isfield(Options.Solver.MCMC,'Proposal')
                if isfield(Options.Solver.MCMC.Proposal,'propDist')
                    error('The prior needs to be specified, if uq_input object is used in sampler proposal')
                end
            end
        end
    end
    % get input from uqlab session
    Opt.Value = uq_getInput;
elseif ~isa(Opt.Value,'uq_input')
  error('The prior is invalid');
end

Internal.Prior = Opt.Value;
Internal.nModelParams = length(Internal.Prior.Marginals);

if ~CUSTOM_LIKELI
    %% MODEL
    [Opt,Options] = uq_process_option(Options,'ForwardModel');
    % missing & invalid
    if Opt.Missing
        % Retrieve model from UQLab and assign
        currModel = uq_getModel;
        if ~isempty(currModel)
            Opt.Value.Model = currModel;
        else
            error('The model is missing');
        end
    else
        if isa(Opt.Value, 'uq_model')
            % single model passed in ForwardModel field
            ForwardModel = Opt.Value; Opt = rmfield(Opt,'Value');
            Opt.Value.Model = ForwardModel;
        else
            % model(s) passed in ForwardModel struct array
            for ii = 1:length(Opt.Value)
                if ~isa(Opt.Value(ii).Model,'uq_model')
                    % supplied forward model is not a uq_model
                    error('The supplied model is not a uq_model');
                end 
            end
        end
    end
    
    % Check PMap
    if length(Opt.Value) == 1
        % only single model supplied - use default parameter map
        Opt.Value.PMap = 1:Internal.nModelParams;
    else
        % check if PMap is supplied
        if isfield(Opt.Value,'PMap')
            % PMap supplied, check if it is consistent with input
            PMapCombined = [Opt.Value.PMap];
            if unique(PMapCombined) ~= Internal.nModelParams
                error('Provided PMap is not consistent with supplied model prior')
            end
        else
            % No PMap supplied, create one by assuming every model takes 
            % the same input
            for ii = 1:length(Opt.Value)
                Opt.Value(ii).PMap = 1:Internal.nModelParams;
            end
        end
    end
    
    Internal.ForwardModel = Opt.Value;
    Internal.nForwardModels = length(Opt.Value);
end

%% DATA
[Opt,Options] = uq_process_option(Options,'Data');
% missing & invalid
if Opt.Missing
  error('The data are missing');
end

if ~isfield(Opt.Value,'y') %data are missing
    error('No data supplied in data group')
end

if ~CUSTOM_LIKELI
    % loop over data groups and perform consistency checks
    % use running index for MOMap, so that by default consecutive model 
    % outputs are adressed by the supplied data
    runI = 0;
    for ii = 1:length(Opt.Value)
        [~,NoutCurr] = size(Opt.Value(ii).y);
        if ~isfield(Opt.Value(ii),'MOMap') % no MOMap - generate
            if length(Internal.ForwardModel) > 1
                % can't handle multiple models without MOMap
                error('In case of multiple models, the MOMap has to be specified fully')
            else
                %generate MOMap for single model (number consecutively)
                Opt.Value(ii).MOMap = [ones(1,NoutCurr);runI+(1:NoutCurr)];
                runI = runI+NoutCurr; %update running index
            end
        elseif size(Opt.Value(ii).MOMap,1) == 1
            % If only one row is given, this refers to the model output id
            if length(Internal.ForwardModel) > 1
                % can't handle multiple models with one line MOMap
                error('In case of multiple models, the MOMap has to be specified fully')
            else
                % put MOMap for single model in matrix
                Opt.Value(ii).MOMap = [ones(1,NoutCurr);Opt.Value(ii).MOMap];
            end
        else
            % check size
            if NoutCurr ~= size(Opt.Value(ii).MOMap,2)
                error('MOMap is not consistent with supplied data')
            end
        end
    end
    
    % check that every output has one unique discrepancy assigned to it
    CombinedMOMap = [Opt.Value.MOMap];
    for ii = 1:Internal.nForwardModels
        CurrModelOMap = CombinedMOMap(2,CombinedMOMap(1,:)==ii);
        if length(unique(CurrModelOMap)) ~= length(CurrModelOMap)
            error('The supplied MOMap does not uniquely address every model output')
        end
    end
else
    if length(Opt.Value) > 1
        error('Multiple discrepancy groups not supported for custom likelihood')
    end
end
Internal.nDataGroups = length(Opt.Value);
Internal.Data = Opt.Value;
uq_addprop(CurrentAnalysis,'Data',Opt.Value);


%% DISCREPANCY DEFAULTS
% Default options for discrepancy
% use unknown sigma with uniform prior from 0 to the square of the mean of
% the observations for each 
if ~CUSTOM_LIKELI
    for ii = 1:Internal.nDataGroups
        DataPointsCurr = Internal.Data(ii).y;
        % clear DiscrepancyPriorOpt and define new
        clear DiscrepancyPriorOpt
        DiscrepancyPriorOpt.Name = sprintf('Prior of sigma %i',ii);
        % for each data dimension an individual discrepancy variance is
        % inferred
        for jj = 1:size(DataPointsCurr,2)
            DiscrepancyPriorOpt.Marginals(jj).Name = 'Sigma2';
            DiscrepancyPriorOpt.Marginals(jj).Type = 'Uniform';
            DiscrepancyPriorOpt.Marginals(jj).Parameters = [0 mean(DataPointsCurr(:,jj)).^2];
        end
        defaultDiscrepancyPrior = uq_createInput(DiscrepancyPriorOpt,'-private');
        % create default discrepancy structure
        DiscrepancyDefault(ii).Type = 'Gaussian';
        DiscrepancyDefault(ii).Prior = defaultDiscrepancyPrior;
        DiscrepancyDefault(ii).Parameters = [];
        % if a parameter is given, in current data group, set prior default to
        % []
        if isfield(Options,'Discrepancy')
            if isfield(Options.Discrepancy(ii),'Parameters')
                if ~isempty(Options.Discrepancy(ii).Parameters)
                    DiscrepancyDefault(ii).Prior = [];
                end
            end
        end
    end
end

%% DISCREPANCY
if ~CUSTOM_LIKELI
    [Opt,Options] = uq_process_option(Options,'Discrepancy',DiscrepancyDefault);
    % missing & invalid
    if and(Opt.Missing, Internal.Display > 0)
        fprintf('The discrepancy was not specified,\nusing unknown i.i.d. Gaussian discrepancy...\n');
    end
    
    % Assign Discrepancy to Internal
    Internal.Discrepancy = Opt.Value;
    
    % Go over data groups and make sure the discrepancies are specified for each data
    % point
    for ii = 1:Internal.nDataGroups
        currDiscrepancyOpt = Opt.Value(ii);
        currDataSize = size(Internal.Data(ii).y);
        
        % get data group size
        switch lower(currDiscrepancyOpt.Type)
          case 'gaussian'
              % covariance matrix
                if isfield(currDiscrepancyOpt,'Parameters')
                    if ~isempty(currDiscrepancyOpt.Parameters)
                        if isfield(currDiscrepancyOpt,'Prior')
                            if ~isempty(currDiscrepancyOpt.Prior)
                                error('Parameters and Prior specified for discrepancy. Only one can be defined.')
                            end
                        end
                        % known discrepancy
                        if isnumeric(currDiscrepancyOpt.Parameters) 
                            ParamSize = size(currDiscrepancyOpt.Parameters);
                            if all(ParamSize==[1,1]) %scalar
                                Internal.Discrepancy(ii).ParamFamily = 'Scalar';
                            elseif and(ParamSize(1)==1,ParamSize(2)==currDataSize(2)) %row
                                Internal.Discrepancy(ii).ParamFamily = 'Row';
                            elseif all(ParamSize==[currDataSize(2),currDataSize(2)]) %matrix
                                %compute cholesky decomposition to make
                                %sure matrix is positive definite
                                [~,p] = chol(currDiscrepancyOpt.Parameters);
                                if and(p == 0,issymmetric(currDiscrepancyOpt.Parameters))
                                    %matrix is positive definite
                                    Internal.Discrepancy(ii).ParamFamily = 'matrix';
                                else
                                    error('Discrepancy matrix is not a covariance matrix')
                                end
                            else
                              error('The discrepancy parameter size is inconsistent');
                            end
                            %update discrepancy info
                            Internal.Discrepancy(ii).ParamType = 'Gaussian';
                            Internal.Discrepancy(ii).nParams = 0;
                            Internal.Discrepancy(ii).ParamKnown = true;
                        else
                            error('Discrepancy parameter has to be numeric')
                        end
                    end
                end
                if isfield(currDiscrepancyOpt,'Prior')
                    if ~isempty(currDiscrepancyOpt.Prior)
                        if isfield(currDiscrepancyOpt,'Parameters')
                            if ~isempty(currDiscrepancyOpt.Parameters)
                                error('Parameters and Prior specified for discrepancy. Only one can be defined.')
                            end
                        end
                        % unknown discrepancy
                        if isa(currDiscrepancyOpt.Prior,'uq_input') 
                            ParamElems = length(currDiscrepancyOpt.Prior.Marginals);
                            %check support of discrepancy prior
                            for jj = 1:ParamElems
                                lowBound = uq_all_invcdf(0, currDiscrepancyOpt.Prior.Marginals(jj));
                                if lowBound < 0
                                    error('Only distributions with positive support can be used as priors for the discrepancy variance')
                                end
                            end
                            if ParamElems==1 %scalar
                                Internal.Discrepancy(ii).ParamFamily = 'Scalar';
                            elseif ParamElems==currDataSize(2) %row
                                Internal.Discrepancy(ii).ParamFamily = 'Row';
                            else
                              error('Discrepancy parameter type is not supported');
                            end
                            % assign discrepancy to internal
                            Internal.Discrepancy(ii).ParamType = 'Gaussian';
                            Internal.Discrepancy(ii).nParams = numel(currDiscrepancyOpt.Prior.nonConst);
                            Internal.Discrepancy(ii).ParamKnown = false;
                        else
                            error('Discrepancy prior has to be UQLab INPUT object')
                        end
                    end
                else
                    error('Either parameters or a prior have to be defined for discrepancy options')
                end
          otherwise
                error('Non supported discrepancy type specified');
        end
    end
    % Check that every data group has an associated discrepancy group
    if ~isequal(length(Internal.Data),length(Internal.Discrepancy))
        error('Number of discrepancy and data groups non consistent')
    end
else
    %custom likelihood
    Internal.Discrepancy(1).nParams = 0;
    Internal.Discrepancy(1).ParamKnown = true;
end
%% FULL PRIOR (including discrepancy parameters and excluding constants)
% copy the prior of the model parameters
PriorOpt.Name = 'Full prior distribution';
PriorOpt.Marginals = rmfield(Internal.Prior.Marginals,'Moments');
PriorOpt.Copula = Internal.Prior.Copula;

% indices of the non constant variables
idNonConst = Internal.Prior.nonConst;	
nTotalParams = Internal.nModelParams; % all parameters counter

% paramDiscrepancyId connects parameters with model or discrepancy parameters 
% 0...model parameters
% i...ith discrepancy parameter
paramDiscrepancyId = zeros(1,nTotalParams); 

% append the discrepancy distribution to the prior
% loop over data groups
for ii = 1:Internal.nDataGroups
    discrepancyCurr = Internal.Discrepancy(ii);
    % number of discrepancy elements
    nDiscrepancyParams = discrepancyCurr.nParams;
    % does the current group have an unknown discrepancy
    if ~discrepancyCurr.ParamKnown
        % add the marginal and copula
        for jj = 1:nDiscrepancyParams
            PriorOpt.Marginals(end+1) = orderfields(...
                rmfield(discrepancyCurr.Prior.Marginals(jj), 'Moments'),...
                PriorOpt.Marginals);
            PriorOpt.Copula.Parameters(end+1,end+1) = 1;
            PriorOpt.Copula.RankCorr(end+1,end+1) = 1;
            % add non constant indices
            if ~strcmp(discrepancyCurr.Prior.Marginals(jj).Type,'Constant')
                idNonConst = [idNonConst,(nTotalParams+1)];
            end
            % update total parameter counter
            nTotalParams = nTotalParams +1;
            paramDiscrepancyId = [paramDiscrepancyId, ii];
        end
    end
end

% add paramDiscrepancyId to internal
Internal.paramDiscrepancyID = paramDiscrepancyId;

% indices and values of the constants
sizePrior = length(PriorOpt.Marginals);
idConst = zeros(1,sizePrior-length(idNonConst));
valConst = zeros(1,sizePrior-length(idNonConst));
jj = 1;
for ii = 1:sizePrior
  if ~any(ii==idNonConst)
    % indices of the constants
    idConst(jj) = ii;
    % values of the constants
    valConst(jj) = PriorOpt.Marginals(ii).Parameters;
    % update index
    jj = jj+1;
  end
end
% constant parameters
								
if ~isempty(idConst)
  PriorOpt.Marginals(idConst) = [];
  PriorOpt.Copula.Parameters(idConst,:) = [];
  PriorOpt.Copula.Parameters(:,idConst) = [];
  PriorOpt.Copula.RankCorr(idConst,:) = [];
  PriorOpt.Copula.RankCorr(:,idConst) = [];
end

% create prior uqlab input object
Internal.FullPrior = uq_createInput(PriorOpt,'-private');


% add property
uq_addprop(CurrentAnalysis,'PriorDist',Internal.FullPrior);
uq_addprop(CurrentAnalysis,'Prior',@(x) uq_eval_pdf(x, Internal.FullPrior));
uq_addprop(CurrentAnalysis,'LogPrior',@(x) uq_eval_logpdf(x, Internal.FullPrior));

% add constant info to analysis object
ConstInfo.idConst = idConst;
ConstInfo.idNonConst = idNonConst;
ConstInfo.constVal = valConst;
Internal.ConstInfo = ConstInfo;

%% AUGMENT INPUT (CONSTANTS)
if isempty(idConst) % without constants
  Aug = @(x) x;
else % with constants
  Aug = @(x) uq_inversion_AugmentInput(x,idNonConst,idConst,valConst);
end

if ~CUSTOM_LIKELI
    %% FORWARD MODEL
    if isempty(idConst) % without constants
      ForwardModel = Internal.ForwardModel;
      %just copy no constant version
      Internal.ForwardModel_WithoutConst = ForwardModel;
    else % with constants, remove constants
      %loop over forward models
      for ii = 1:Internal.nForwardModels
        ModelOpt.Name = Internal.ForwardModel(ii).Model.Name;
        ModelOpt.isVectorized = true;
        ModelOpt.mHandle = @(x) uq_evalModel(Internal.ForwardModel(ii).Model,Aug(x));
        ForwardModel(ii).Model = uq_createModel(ModelOpt,'-private');
        %update PMap
        PMap = Internal.ForwardModel(ii).PMap;
        idConstCurr = ConstInfo.idConst;
        for jj = 1:length(idConstCurr)
            %remove const id
            PMap = PMap(PMap~=idConstCurr(jj));
            %reduce index for larger id
            PMap(PMap > idConstCurr(jj)) = PMap(PMap > idConstCurr(jj)) - 1;
            idConstCurr(jj+1:end) = idConstCurr(jj+1:end) - 1;
        end
        ForwardModel(ii).PMap = PMap;
        %assign to internal
        Internal.ForwardModel_WithoutConst = ForwardModel;
      end
    end
    uq_addprop(CurrentAnalysis,'ForwardModel',ForwardModel);
end

%% LIKELIHOOD
% Loop over individual data groups and specify each groups likelihood
% functions

if ~CUSTOM_LIKELI
    % Loop over data groups
    for ii = 1:Internal.nDataGroups
        discrepancyCurr = Internal.Discrepancy(ii);
        dataPointsCurr = Internal.Data(ii).y; % measurements
        % switch discrepancy type
        switch lower(discrepancyCurr.ParamType) 
            case 'gaussian'
                switch lower(discrepancyCurr.ParamFamily) % switch discrepancy family
                    case 'scalar'
                        Likelihood_g = @(modelRuns,param) ...
                            uq_inversion_normpdf(modelRuns,dataPointsCurr,param,'Scalar');
                        LogLikelihood_g = @(modelRuns,param) ...
                            uq_inversion_lognormpdf(modelRuns,dataPointsCurr,param,'Scalar');
                        LikelihoodSamples_g = @(y,param) normrnd(y,sqrt(param));
                    case 'row'
                        Likelihood_g = @(modelRuns,param) ...
                            uq_inversion_normpdf(modelRuns,dataPointsCurr,param,'Row');
                        LogLikelihood_g = @(modelRuns,param) ...
                            uq_inversion_lognormpdf(modelRuns,dataPointsCurr,param,'Row');
                        LikelihoodSamples_g = @(y,param) normrnd(y,sqrt(param));
                    case 'matrix'
                        Likelihood_g = @(modelRuns,param) ...
                            uq_inversion_normpdf(modelRuns,dataPointsCurr,param,'Matrix');
                        LogLikelihood_g = @(modelRuns,param) ...
                            uq_inversion_lognormpdf(modelRuns,dataPointsCurr,param,'Matrix');
                        LikelihoodSamples_g = @(y,param) mvnrnd(y,param);
                    otherwise
                        error('The discrepancy parameter has the wrong size');
                end
            otherwise
                error('The discrepancy type is not supported');
        end
        %add properties
        Internal.Discrepancy(ii).likelihood_g = Likelihood_g;
        Internal.Discrepancy(ii).logLikelihood_g = LogLikelihood_g;
        Internal.Discrepancy(ii).likelihoodSamples_g = LikelihoodSamples_g;
    end
    %now add a global likelihood evaluator
    uq_addprop(CurrentAnalysis,'Discrepancy',Internal.Discrepancy);

    Internal.Likelihood = @(x) uq_inversion_likelihood(Aug(x),Internal,'Likelihood');
    Internal.LogLikelihood = @(x) uq_inversion_likelihood(Aug(x),Internal,'LogLikelihood');

    uq_addprop(CurrentAnalysis,'Likelihood',Internal.Likelihood);
    uq_addprop(CurrentAnalysis,'LogLikelihood',Internal.LogLikelihood);
    
    %add posterior handle
    uq_addprop(CurrentAnalysis,'UnnormPosterior',@(x) Internal.Likelihood(x).*uq_eval_pdf(x, Internal.FullPrior));
    uq_addprop(CurrentAnalysis,'UnnormLogPosterior',@(x) Internal.LogLikelihood(x) + uq_eval_logpdf(x, Internal.FullPrior));
else
    %custom likelihood
    Internal.LogLikelihood = @(x) Internal.LogLikelihood(Aug(x),Internal.Data.y);
    
    %add field
    uq_addprop(CurrentAnalysis,'LogLikelihood',Internal.LogLikelihood);
    
    %add posterior handle
    uq_addprop(CurrentAnalysis,'UnnormLogPosterior',@(x) Internal.LogLikelihood(x) + uq_eval_logpdf(x, Internal.FullPrior));
end

%% SOLVER DEFAULTS
% Default options for MCMC 
SolverDefault.Type = 'MCMC';
SolverDefault.MCMC.Sampler = 'AIES';
SolverDefault.MCMC.StoreModel = true;
SolverDefault.MCMC.Visualize.Parameters = 0;
SolverDefault.MCMC.Visualize.Interval = 0;
% Default options per sampler
% MH
% 3000 iterations, 10 Chains, 10% of prior covariance
SolverDefaultsSamplerMH.Steps = 3000;
SolverDefaultsSamplerMH.NChains = 10;
SolverDefaultsSamplerMH.Proposal.PriorScale = 0.1;
% AM
% 3000 iterations, 300 with initial covariance matrix, 10 Chains,
SolverDefaultsSamplerAM.Steps = 3000;
SolverDefaultsSamplerAM.NChains = 10;
SolverDefaultsSamplerAM.Proposal.PriorScale = 0.1;
SolverDefaultsSamplerAM.T0 = 300;
SolverDefaultsSamplerAM.Epsilon = 1e-6;
% HMC
% 500 iterations, 10 Chains, mass proportional to second moment, 
% 50 leapfrog steps, leapfrog size 
SolverDefaultsSamplerHMC.Steps = 300;
SolverDefaultsSamplerHMC.NChains = 10;
SolverDefaultsSamplerHMC.Mass = 1;
SolverDefaultsSamplerHMC.LeapfrogSteps = 10;
SolverDefaultsSamplerHMC.LeapfrogSize = 0.01;
% AIES
% 300 iterations, 100 chains, a = 2
SolverDefaultsSamplerAIES.Steps = 300;
SolverDefaultsSamplerAIES.NChains = 100;
SolverDefaultsSamplerAIES.a = 2;

%% SOLVER
% First level defaults
[Opt,Options] = uq_process_option(Options,'Solver',SolverDefault);
% missing & invalid
if Opt.Invalid
  error('The solver is invalid');
elseif and(Opt.Missing, Internal.Display > 0)
  fprintf('The solver was not specified, using MCMC');
end

% Second level defaults for MCMC Sampler
if strcmp(Opt.Value.Type,'MCMC')
    % Sampler
    SamplerOpt = uq_process_option(Opt.Value.MCMC,'Sampler',SolverDefault.MCMC.Sampler);
    % missing & invalid
    if SamplerOpt.Invalid
      error('The sampler is invalid');
    elseif and(SamplerOpt.Missing, Internal.Display > 0)
      fprintf('The sampler was not specified, using affine invariant ensemble sampler');
    end
    % assign second level to first level
    Opt.Value.MCMC.Sampler = SamplerOpt.Value;
    
    % StoreModel
    StoreModelOpt = uq_process_option(Opt.Value.MCMC,'StoreModel',SolverDefault.MCMC.StoreModel);
    % missing & invalid
    if StoreModelOpt.Invalid
      error('The model storage options are invalid');
    end
    % assign second level to first level
    Opt.Value.MCMC.StoreModel = StoreModelOpt.Value;
    % If custom likelihood is supplied, turn off storage of model
    % evaluations
    if CUSTOM_LIKELI
        Opt.Value.MCMC.StoreModel = false;
    end
    
    % Visualize
    VisualizeOpt = uq_process_option(Opt.Value.MCMC,'Visualize',SolverDefault.MCMC.Visualize);
    % missing & invalid
    if VisualizeOpt.Invalid
      error('The visualization options are invalid');
    end
    % assign second level to first level
    Opt.Value.MCMC.Visualize = VisualizeOpt.Value;
end


%run through solvers and check if all options are set
switch upper(Opt.Value.Type)
    case 'NONE'
        %add likelihood to solver field
        Opt.Value.LogPrior = @(x) uq_eval_logpdf(x,Internal.FullPrior);
        Opt.Value.LogLikelihood = Internal.LogLikelihood;
    case 'MCMC'
        % are all general MCMC options set
        if ~isfield(Opt.Value.MCMC,'Sampler')
            error('The sampler was not specified');
        end
        %switch between samplers
        switch upper(Opt.Value.MCMC.Sampler)
            %are all sampler specific options set
            case 'MH'
                %defaults
                MCMCOpt = uq_process_option(Opt.Value,'MCMC',SolverDefaultsSamplerMH);
                Opt.Value.MCMC = MCMCOpt.Value;
                %metropolis hastings
                if ~isfield(Opt.Value.MCMC,'Steps')
                    error('Did not set the number of steps');
                elseif ~isfield(Opt.Value.MCMC,'Proposal')
                    error('Did not set proposal properties');
                end
                %Initialize proposal distribution
                PriorMoments = reshape([Internal.FullPrior.Marginals.Moments],2,[]);
                PriorVariance = PriorMoments(2,:).^2;
                Opt.Value.MCMC.Proposal = uq_initialize_uq_inversion_proposal(Opt.Value.MCMC.Proposal,PriorVariance);
            case 'AM'
                %defaults
                MCMCOpt = uq_process_option(Opt.Value,'MCMC',SolverDefaultsSamplerAM);
                Opt.Value.MCMC = MCMCOpt.Value;
                %adaptive metropolis
                if ~isfield(Opt.Value.MCMC,'Steps')
                    error('Did not set the number of steps');
                elseif ~isfield(Opt.Value.MCMC,'T0')
                    error('Did not set a stepcount for an initial covariance');
                elseif ~isfield(Opt.Value.MCMC,'Proposal')
                    error('Did not set proposal properties');
                elseif ~isfield(Opt.Value.MCMC,'Epsilon')
                    error('Did not specify epsilon value')
                end
                if length(Opt.Value.MCMC.Epsilon) == numel(Opt.Value.MCMC.Epsilon)
                    if ~and(size(Opt.Value.MCMC.Epsilon,1) == 1,size(Opt.Value.MCMC.Epsilon,2) == 1)
                        %is not scalar
                        error('Epsilon has to be a scalar')
                    end
                else
                    error('Epsilon has to be given as a scalar or a vector')
                end
                %Initialize proposal distribution
                PriorMoments = reshape([Internal.FullPrior.Marginals.Moments],2,[]);
                PriorVariance = PriorMoments(2,:).^2;
                Opt.Value.MCMC.Proposal = uq_initialize_uq_inversion_proposal(Opt.Value.MCMC.Proposal,PriorVariance);
            case 'HMC'
                %defaults
                MCMCOpt = uq_process_option(Opt.Value,'MCMC',SolverDefaultsSamplerHMC);
                Opt.Value.MCMC = MCMCOpt.Value;
                %hamiltonian monte carlo
                if ~isfield(Opt.Value.MCMC,'Steps')
                    error('Did not set the number of steps');
                elseif ~isfield(Opt.Value.MCMC,'LeapfrogSteps')
                    error('Did not set a number for leap frog steps');
                elseif ~isfield(Opt.Value.MCMC,'LeapfrogSize')
                    error('Did not set a step size specifier');
                elseif ~isfield(Opt.Value.MCMC,'Mass')
                    error('Did not set the mass');
                end
                %distinguish between scalar and matrix case
                mass = Opt.Value.MCMC.Mass;
                if and(all(size(mass)==1), mass > 0)
                    %positive scalar
                    Opt.Value.MCMC.Mass = mass*eye(length(ConstInfo.idNonConst));
                else 
                    %compute cholesky decomposition to make
                    %sure mass matrix is positive definite
                    [~,p] = chol(mass);
                    if ~and(p == 0,issymmetric(Opt.Value.MCMC.Mass))
                        error('Supplied mass matrix is not positive definite')
                    end
                end
            case 'AIES'
                %defaults
                MCMCOpt = uq_process_option(Opt.Value,'MCMC',SolverDefaultsSamplerAIES);
                Opt.Value.MCMC = MCMCOpt.Value;
                %affine invariant ensemble sampler
                if ~isfield(Opt.Value.MCMC,'Steps')
                    error('Did not set the number of steps');
                elseif ~isfield(Opt.Value.MCMC,'a')
                    error('Did not set the a parameter');
                end
            otherwise
                error('The specified sampler is not supported');
        end
        
        % draw NChains number of Seeds from prior distribution
        if ~isfield(Opt.Value.MCMC,'Seed')
            %No initial point is set - check if number of chains is given
            if ~isfield(Opt.Value.MCMC,'NChains')
                error('Did not specify intial points or number of chains for MCMC sampler')
            else
                %number of chains set - draw appropriate number of samples
                %from prior
                Seed = uq_getSample(Internal.FullPrior,Opt.Value.MCMC.NChains,'Sobol')';
                Opt.Value.MCMC.Seed = Seed;
            end
        end
        
        %if seeds are given in 2d-array restructure to 3d-array
        if ismatrix(Opt.Value.MCMC.Seed)
            Seed = Opt.Value.MCMC.Seed;
            Opt.Value.MCMC.Seed = reshape(Seed,1,size(Seed,1),size(Seed,2));
        end
        
        if ~(size(Opt.Value.MCMC.Seed,2)==length(Internal.FullPrior.Marginals))
            error('MCMC seeds do not match prior distribution size')
        end
        
        %taking care of constants in seeds
        if size(Opt.Value.MCMC.Seed,2) > numel(ConstInfo.idNonConst)
            if Internal.Display > 0
                fprintf('Constants given in MCMC seed, removing constants...')
            end
            Opt.Value.MCMC.Seed = Opt.Value.MCMC.Seed(:,ConstInfo.idNonConst,:);
        end
        
        %add likelihood to solver field
        Opt.Value.LogPrior = @(x) uq_eval_logpdf(x,Internal.FullPrior);
        Opt.Value.LogLikelihood = Internal.LogLikelihood;
    otherwise
        error('The specified solver is not supported');
end

Internal.Solver = Opt.Value;

%% UNPROCESSED OPTIONS
if isfield(Options,'Type')
  Options = rmfield(Options,'Type');
end
if isfield(Options,'Name')
  Options = rmfield(Options,'Name');
end
uq_options_remainder(Options);

%% ACTUAL OPTIONS
CurrentAnalysis.Internal = Internal;

%% RUN ANALYSIS
uq_runAnalysis(CurrentAnalysis);

%% SUCCESS
Success = 1;
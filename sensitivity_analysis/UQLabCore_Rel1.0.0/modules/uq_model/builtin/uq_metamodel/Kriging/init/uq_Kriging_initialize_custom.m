function success = uq_Kriging_initialize_custom( current_model, Options )
% success = UQ_KRIGING_INITIALIZE_CUSTOM(CURRENT_MODEL,OPTIONS) initialize a
%     custom Kriging model (predictor only) in CURRENT_MODEL based on the options given
%     in the OPTIONS structure
%
% See also: UQ_KRIGING_INITIALIZE

success = 0;

assert(isfield(Options.Kriging,'beta') && ...
    ~any(isempty([Options.Kriging(:).beta])),...
    'For Custom Kriging the beta field is needed!');

assert(isfield(Options.Kriging,'Trend') && ...
    ~any(isempty([Options.Kriging(:).Trend])),...
    'For Custom Kriging the Trend field is needed!');

% make sure that also the field Trend.Type has been set
for ii = 1 : length(Options.Kriging)
    assert(all(isfield(Options.Kriging(ii).Trend,'Type')) && ...
        ~isempty([Options.Kriging(ii).Trend.Type]),...
        'The Trend.Type field is needed!');
end
assert(isfield(Options.Kriging,'sigmaSQ') &&...
    ~any(isempty([Options.Kriging(:).sigmaSQ])),...
    'For Custom Kriging the sigmaSQ field is needed!');

assert(isfield(Options.Kriging,'theta') && ...
    ~any(isempty([Options.Kriging(:).theta])),...
    'For Custom Kriging the theta field is needed!');

assert(isfield(Options.Kriging,'Corr') && ...
    ~any(isempty([Options.Kriging(:).Corr])),...
    'For Custom Kriging the Corr(correlation function options)field is needed!');

Nout = length(Options.Kriging);
NuggetValue = cell(Nout,1);

for ii = 1 : Nout
    assert(all(isfield(Options.Kriging(ii).Corr,'Family')) && ...
        ~isempty([Options.Kriging(ii).Corr.Family]),...
        'The correlation family field (.Corr.Family) is needed!');
    assert(all(isfield(Options.Kriging(ii).Corr,'Type')) && ...
        ~isempty([Options.Kriging(ii).Corr.Type]),...
        'The correlation type field (.Corr.Type) is needed!');
    % the Isotropic flag is optional
    if isfield(Options.Kriging(ii).Corr,'Isotropic') && ...
            ~isempty([Options.Kriging(ii).Corr.Isotropic])
        corrIsIsotropic(ii) = Options.Kriging(ii).Corr.Isotropic ;
    else
        % by default the correlation function is considered
        % anisotropic
        corrIsIsotropic(ii) = false ;
    end
    
    % get the Nugget or set its default value
    if isfield(Options.Kriging(ii).Corr,'Nugget') && ...
            ~isempty([Options.Kriging(ii).Corr.Nugget])
        NuggetValue{ii} = Options.Kriging(ii).Corr.Nugget ;
    else
        % by default the a small nugget value is set
        % (the same default value is used during typical Kriging
        % initialization)
        NuggetValue{ii} = 1e-10 ;
    end
end


assert(isfield(Options,'ExpDesign'),...
    'For Custom Kriging the Experimental Design X and Y need to be defined!');


assert(isfield(Options.ExpDesign,'X') && ~isempty([Options.ExpDesign.X]),...
    'For Custom Kriging the Experimental Design X field is needed!');


assert(isfield(Options.ExpDesign,'Y') && ~isempty([Options.ExpDesign.Y]),...
    'For Custom Kriging the Experimental Design Y field is needed!');

M = size(Options.ExpDesign.X, 2) ;
% M = length(Options.Input.nonConst);
current_model.Internal.Runtime.M = M;
try
    current_model.Internal.Runtime.nonConstIdx = current_model.Options.Input.nonConst;
    current_model.Internal.Runtime.MnonConst = length(current_model.Options.Input.nonConst);
catch
    current_model.Internal.Runtime.nonConstIdx = 1:M;
    current_model.Internal.Runtime.MnonConst  = M;
end

% Make sure that the number of outputs is consistent in options and
% ExpDesign.Y
assert(isequal(size(Options.ExpDesign.Y,2), length(Options.Kriging)) , ...
    'The length of Options.Kriging is not equal to the number of outputs as defined in the Experimental Design!');

% Make sure that X and Y matrices have equal length (i.e. refer to
% the same number of samples)

assert (isequal(size(Options.ExpDesign.X,1), size(Options.ExpDesign.Y,1)), ...
    'The length of X and Y of the Experimental Design is not equal!')

% Get the number of outputs
Nout = length(Options.Kriging);


% assign an empty input module
current_model.Internal.Input = [];

% now start building the custom metamodel
uq_addprop(current_model, 'ExpDesign');
current_model.ExpDesign.X = Options.ExpDesign.X ;
current_model.ExpDesign.U = Options.ExpDesign.X ;

% Filter out only the degrees of freedom corresponding to non-constants:
nonConstIdx = current_model.Internal.Runtime.nonConstIdx;
U = current_model.ExpDesign.U(:,nonConstIdx);

current_model.ExpDesign.Y = Options.ExpDesign.Y ;
current_model.Internal.ExpDesign.varY = var(current_model.ExpDesign.Y);
current_model.ExpDesign.NSamples = size(U, 1);

current_model.Internal.Scaling = 0;
current_model.Internal.Runtime.Nout = Nout;

for oo = 1 : Nout
    % runtime fields
    current_model.Internal.Runtime.current_output = oo ;
    
    % trend fields
    
    % if no truncation options have been set do not truncate
    % anything from the polynomial basis
    if ~isfield(Options.Kriging(oo).Trend, 'TruncOptions') || ...
            isempty(Options.Kriging(oo).Trend)
        Options.Kriging(oo).Trend.TruncOptions.qNorm = 1 ;
    end
    % if no polytypes have been set assume simple polynomial basis
    if ~isfield(Options.Kriging(oo).Trend, 'PolyTypes') || ...
            isempty(Options.Kriging(oo).Trend.PolyTypes)
        Options.Kriging(oo).Trend.PolyTypes = repmat({'simple_poly'},M,1) ;
    end
    
    % Trend
    current_model.Internal.Kriging(oo).Trend = ...
        uq_Kriging_initializeTrend( Options.Kriging(oo).Trend, M);
    % for now the default trend handle is going to be
    %  assigned always, in case of custom kriging:
    current_model.Internal.Kriging(oo).Trend.Handle = @uq_Kriging_eval_F;
    current_model.Internal.Kriging(oo).Trend.beta = ...
        Options.Kriging(oo).beta;
    current_model.Internal.Kriging(oo).Trend.F = ...
        uq_Kriging_eval_F( U, current_model );
    
    % GP fields
    current_model.Internal.Kriging(oo).GP.Corr.Type = ...
        Options.Kriging(oo).Corr.Type;
    current_model.Internal.Kriging(oo).GP.Corr.Family = ...
        Options.Kriging(oo).Corr.Family;
    current_model.Internal.Kriging(oo).GP.Corr.Isotropic = ...
        corrIsIsotropic(oo);
    % for now the default covariance handle is going to be
    %      assigned always, in case of custom kriging:
    current_model.Internal.Kriging(oo).GP.Corr.Handle = ...
        @uq_eval_Kernel;
    current_model.Internal.Kriging(oo).GP.Corr.Nugget = NuggetValue{oo} ;
    current_model.Internal.Kriging(oo).GP.sigmaSQ = ...
        Options.Kriging(oo).sigmaSQ;
    current_model.Internal.Kriging(oo).GP.R = ...
        uq_eval_Kernel(U,U, ...
        Options.Kriging(oo).theta, ...
        current_model.Internal.Kriging(oo).GP.Corr ) ;
    
    % optim fields
    current_model.Internal.Kriging(oo).Optim.Theta = ...
        Options.Kriging(oo).theta ;
    
    % set an empty cached field
    current_model.Internal.Kriging(oo).Cached = [];
end


% By default enable the KeepCache feature so that subsequent
% predictions are faster (that is from the 2nd call onwards)
current_model.Internal.KeepCache = 1;
% add the remaining runtime arguments that may be needed
current_model.Internal.Runtime.isCalculated = true;

% done
success = 1;




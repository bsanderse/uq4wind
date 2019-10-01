function varargout = uq_eval_uq_default_model(current_model,X)
% UQ_EVAL_UQ_DEFAULT_MODEL evaluates a model defined in an m-file, m-string
% or m-handle
%
% See also: UQ_EVAL_UQ_METAMODEL, UQ_INITIALIZE_UQ_DEFAULT_MODEL

% assume 1 output argument when none is selected
num_of_out_args = max(nargout,1);

% do nothing if X is empty
if isempty(X)
    [varargout{1:num_of_out_args}] = deal([]);
    return;
end

% check if parameters exist
PARAMS_EXIST = ~isempty(current_model.Parameters);
% get the parameters if any
if PARAMS_EXIST
    params = current_model.Parameters ;
end

% retrieve the m-file handle
if ~isfield(current_model.Internal, 'fHandle') || ...
        isempty(current_model.Internal.fHandle)
    error('The model does not seem to be properly initialized! (Internal field fHandle is missing.)')
end

model_handle = current_model.Internal.fHandle ;

%% Calculate the model response

outFull=cell(1,num_of_out_args);
outCurr=cell(1,num_of_out_args);
if ~current_model.isVectorized
    % The current model is NOT vectorized
    
    if PARAMS_EXIST
        % there are parameters
        for ii = 1 : size(X,1)
            [outCurr{1:num_of_out_args}]  = model_handle(X(ii,:), params);
            outFull = cellfun(@(x1,x2)cat(1,x1,x2),outFull,outCurr,...
                'UniformOutput',0);
        end
    else
        %there are no parameters
        for ii = 1 : size(X,1)
            [outCurr{1:num_of_out_args}]  = model_handle(X(ii,:));
            outFull = cellfun(@(x1,x2)cat(1,x1,x2),outFull,outCurr,...
                'UniformOutput',0);
        end
    end
    varargout = outFull;
    
else 
    % the mfile is vectorized
    if PARAMS_EXIST
        % there are parameters
        [varargout{1:num_of_out_args}]  = model_handle(X, params);
    else
        %there are no parameters
        [varargout{1:num_of_out_args}]  = model_handle(X);
    end
end
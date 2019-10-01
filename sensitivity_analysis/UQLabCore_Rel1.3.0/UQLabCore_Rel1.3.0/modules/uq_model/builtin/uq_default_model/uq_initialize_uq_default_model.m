function success = uq_initialize_uq_default_model(module)
% UQ_INITIALIZE_UQ_DEFAULT_MODEL initializes the default model in UQLab
%
% See also: UQ_EVAL_UQ_DEFAULT_MODEL


success  = 0;

%% retrieve the current model
if exist('module', 'var')
    current_model = uq_getModel(module);
else
    current_model = uq_getModel;
end


%% RETRIEVE THE OPTIONS AND PARSE THEM
Options = current_model.Options;

% Make sure that some kind of model source is given
isMFILE = isfield(Options, 'mFile') && ~isempty(Options.mFile) ;
isMSTRING = isfield(Options, 'mString') && ~isempty(Options.mString) ;
isMHANDLE = isfield(Options, 'mHandle') && ~isempty(Options.mHandle) ;

if isMFILE + isMSTRING + isMHANDLE > 1 
    error('Multiple model definitions found!');
end

if isMFILE + isMSTRING + isMHANDLE < 1 
    error('The model property mFile, or mString or mHandle needs to be defined!')
end


if isMFILE
    uq_addprop(current_model,'mFile',Options.mFile);
    % Add the function handle field 
    current_model.Internal.fHandle = str2func(Options.mFile) ;
end

if isMSTRING
    % check whether the @(X) exists in the beginning of the string
    % if not add it
    if isempty(strfind(Options.mString,'@(X)')) && ...
            isempty(strfind(Options.mString,'@'))
        uq_addprop(current_model,'mString',['@(X, P)', Options.mString]);
    else
        uq_addprop(current_model,'mString',Options.mString);
    end
    % Add the function handle field 
    current_model.Internal.fHandle = str2func(current_model.mString) ;
end

if isMHANDLE
    uq_addprop(current_model,'mHandle',Options.mHandle);
    % Add the function handle field 
    current_model.Internal.fHandle = Options.mHandle;
end



% Include the parameters if any
if isfield(Options,'Parameters')
    uq_addprop(current_model,'Parameters',Options.Parameters);
else 
    uq_addprop(current_model,'Parameters');
end


% Parse the Vectorize option 
if isMFILE
    defaultIsVectorized = 1;
else
    defaultIsVectorized = 0;
end

[optVect, Options] = uq_process_option(Options, 'isVectorized',...
    defaultIsVectorized, {'logical','double'});
if optVect.Invalid
    warning('The isVectorized option was invalid. Using the default value instead.')
end
if optVect.Missing
    %do nothing, just silently assign the default value
end

uq_addprop(current_model,'isVectorized',optVect.Value);




% Store the mfile location in Internal
if isMFILE
    try
        current_model.Internal.Location = which(Options.mFile);
    catch
        % in case which() gives an error it's due to Options.Function being
        % an anonymous function
        error('The m-file : %s cannot be found!', Options.mFile);
    end
else
    current_model.Internal.Location = [] ;
end

%% Set the isInitialized flag
current_model.Internal.Runtime.isInitialized = true;

success = 1;
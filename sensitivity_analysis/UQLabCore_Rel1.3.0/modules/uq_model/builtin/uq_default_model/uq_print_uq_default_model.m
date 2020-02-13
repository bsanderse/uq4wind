function uq_print_uq_default_model( module, varargin )
%UQ_PRINT_UQ_DEFAULT_MODEL prints out information regarding a model module
% 
% See also: UQ_EVAL_UQ_DEFAULT_MODEL


%% CONSISTENCY CHECKS
if ~strcmp(module.Type, 'uq_default_model')
    error('uq_print_uq_default_model only operates on objects of type ''uq_default_model''')
end

%% PRINT
fprintf('-------------------------------------------\n')
fprintf('Model object name:\t%s\n', module.Name)

if isprop(module, 'mFile') && ~isempty(module.mFile) 
    fprintf('Defined by: \t\t%s (%s)\n', module.mFile,'m-file')
    fprintf('m-file location: \t%s\n', module.Internal.Location)
elseif isprop(module, 'mString') && ~isempty(module.mString) ;
    fprintf('Defined by: \t\t%s (%s)\n', module.mString,'m-string')
elseif isprop(module, 'mHandle') && ~isempty(module.mHandle) ;
    fprintf('Defined by: \t\t%s (%s)\n', func2str(module.mHandle),'m-handle')
else
    error('The given model object does not seem to be properly configured/initialized!')
end

fprintf('Parameters: \t\t[%s]\n', uq_sprintf_mat(module.Parameters))
if module.isVectorized
    fprintf('Vectorized: \t\ttrue\n')
else
    fprintf('Vectorized: \t\tfalse\n')
end

fprintf('-------------------------------------------\n')
function varargout = uq_eval_uq_uqlink(current_model,X,varargin)
% UQ_EVAL_UQ_LINK evaluates a UQLink model
%
% See also:UQ_INITIALIZE_UQ_UQLINK, UQ_WRAPPER_UQLINK

% assume 1 output argument when none is selected
num_of_out_args = max(nargout,1) ;

% do nothing if X is empty
if isempty(X)
    [varargout{1:num_of_out_args}] = deal([]);
    return;
end

% Evaluate the auxiliary model that will call the third-party software
if nargin > 3
    action = varargin{1} ;
    matfile = varargin{2} ;
    [varargout{1:num_of_out_args}] = uq_wrapper_uqlink(X, current_model.Internal, action, matfile) ;
    
elseif nargin > 2
    action = varargin{1} ;
    [varargout{1:num_of_out_args}] = uq_wrapper_uqlink(X, current_model.Internal, action) ;
    
else
    [varargout{1:num_of_out_args}] = uq_wrapper_uqlink(X, current_model.Internal) ;
    
end

end

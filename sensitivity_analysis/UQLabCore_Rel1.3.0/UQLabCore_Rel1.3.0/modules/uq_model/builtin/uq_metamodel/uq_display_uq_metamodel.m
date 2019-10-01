function varargout = uq_display_uq_metamodel(current_model, varargin)
% varargout = UQ_DISPLAY_UQ_METAMODEL(CURRENT_MODEL,VARARGIN): define the
%     behavior of the uq_display function for uq_metamodel objects
%
% See also UQ_PCE_DISPLAY, UQ_KRIGING_DISPLAY

%% Choose the proper display function depending on the metatype
switch lower(current_model.MetaType)
    case 'pce' 
        [varargout{1:nargout}] = uq_PCE_display(current_model, varargin{:});
    case 'kriging'
        [varargout{1:nargout}] = uq_Kriging_display(current_model, varargin{:});
        
    case 'pck'
        if nargin == 1
            uq_PCK_display(current_model);
        else
            uq_PCK_display(current_model, varargin{:});
        end
    case 'lra'
        if nargin == 1
            uq_LRA_display(current_model);
        else
            uq_LRA_display(current_model, varargin{:});
        end
        
    case 'svr'
        if nargin > 1
            uq_SVR_display(current_model, varargin{:});
        else
            uq_SVR_display(current_model);
        end
        
    case 'svc'
        if nargin > 1
            uq_SVC_display(current_model, varargin{:});
        else
            uq_SVC_display(current_model);
        end
end

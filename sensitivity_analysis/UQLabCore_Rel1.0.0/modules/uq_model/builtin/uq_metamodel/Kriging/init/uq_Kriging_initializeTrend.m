function [ Trend, EVTs ] = uq_Kriging_initializeTrend( Trend, M, Input )
% UQ_KRIGING_INITIALIZETREND : initialization of the Kriging trend
%
% See also UQ_KRIGING_INITIALIZE, UQ_KRIGING_EVAL_F
EVTs = [];

switch lower(Trend.Type)
    case 'simple'
        Trend.CustomF = initialize_customF(Trend.CustomF,M) ;
    case 'ordinary'
        Trend.Degree = 0;
        Trend.CustomF = [];
    case 'linear'
        Trend.Degree = 1;
        Trend.CustomF = [];
    case 'quadratic'
        Trend.Degree = 2;
        Trend.CustomF = [];
    case 'polynomial'
        Trend.CustomF = [];
    case 'custom'
        Trend.CustomF = initialize_customF(Trend.CustomF,M) ;
    otherwise
        error('Unknown trend type!')
end
% Check the truncation options
% 1) Check whether a Custom truncation value has been selected. If yes get the
% trend degree.
if isfield(Trend.TruncOptions, 'Custom')
    Trend.Degree = max(sum(Trend.TruncOptions.Custom,1));
end
% 2) Check the validity of the .qNorm definition
if isfield(Trend.TruncOptions, 'qNorm') && (Trend.TruncOptions.qNorm < 0 || Trend.TruncOptions.qNorm > 1)
    
    EVTs.Type = 'W';
    EVTs.Message = ['Invalid TruncOptions Value! ', ...
        'Setting the default Truncation settings instead.'];
    EVTs.eventID = 'uqlab:metamodel:kriging:initialize:trendtruncation_invalid';
    
    Trend.TruncOptions = TrendDefaults.TruncOptions ;
end
% 3) .MaxInteraction option: no checks are taking place 



% Trend polytypes validation and initialization of type 'auto'
polyClass = class(Trend.PolyTypes);
switch lower(polyClass)
    case 'char'
        % A string will only be accepted only if M == 1
        if M > 1 && ~strcmpi(Trend.PolyTypes,'auto')
            error('PolyTypes should be a cell array of strings, with length M = %i !',M)
        end
        if M == 1 && ~strcmpi(Trend.PolyTypes,'auto')
            Trend.PolyTypes = {Trend.PolyTypes};
        end
        if strcmpi(Trend.PolyTypes,'auto')
            Trend.PolyTypes = ...
                uq_auto_retrieve_poly_types(Input);
        end
    case 'cell'
        dimensionCheckOK = sum(size(Trend.PolyTypes) == [M,1])==2 ||...
            sum( size(Trend.PolyTypes) == [1,M])==2;
        if M>1 && ~dimensionCheckOK
            error('PolyTypes should either have a single or M = %i elements!',M)
        end
        
    otherwise
        error('Polytypes should either be a string or a cell array!')
end





end



function  cF = initialize_customF(CustomF,M)

switch class(CustomF)
    case 'char'
        %transform the string to function handle cell
        cF = {str2func(CustomF)} ;
    case 'cell'
        %transform each string to function handle
        for jj = 1 : length(CustomF)
            switch class(CustomF{jj})
                case 'char'
                    cF{jj} = str2func(CustomF{jj}) ;
                case 'function_handle'
                    %do nothing
                    cF{jj} = CustomF{jj} ;
                otherwise
                    error('Only strings and function handles are accepted when CustomF is a cell array!')
            end
        end
    case 'function_handle'
        %transform the string to function handle cell
        cF = {CustomF} ;
    otherwise
        if isnumeric(CustomF)
            dimensionCheckOK = sum(size(CustomF) == [M,1])==2 ||...
                sum( size(CustomF) == [1,M])==2;
            if ~dimensionCheckOK && M > 1
                error('CustomF dimension mismatch compared to Input dimension!')
            end
        else
            error(['Trend function cannot be defined as '...
                class(CustomF),...
                '! Supported types are numeric, string(containing the',...
                ' function name),cell array of strings and function handle(or cell array of function handles).'])
        end
        cF = CustomF ;
end
% end

end
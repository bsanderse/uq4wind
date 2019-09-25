function uq_Kriging_print(KRGModel, outArray, varargin)
% UQ_KRIGING_PRINT(KRGMODEL,OUTARRAY,VARARGIN): pretty print information on the
%     Kriging object in KRGMODEL for the specified set of output components
%     OUTARRAY (default: OUTARRAY = 1).
%
% See also: UQ_KRIGING_DISPLAY,UQ_PCE_PRINT,UQ_PRINT_UQ_METAMODEL


%% Consistency checks and command line parsing
if ~exist('outArray', 'var') || isempty(outArray)
    outArray = 1;
    if length(KRGModel.Kriging) > 1;
        warning('The selected Kriging metamodel has more than 1 output. Only the 1st output will be printed');
        fprintf('You can specify the outputs you want to be displayed with the syntax:\n')
        fprintf('uq_print(KRGModel, OUTARRAY)\nwhere OUTARRAY is the index of desired outputs, e.g. 1:3 for the first three\n\n')
    end
end
if max(outArray) > length(KRGModel.Kriging)
    error('Requested output range is too large') ;
end

%% parsing the residual command line
% initialization
if nargin > 2
    parse_keys = {'beta', 'theta','F','optim','trend','gp','GP','R'};
    parse_types = {'f', 'f','f', 'f','f', 'f','f','f'};
    [uq_cline, varargin] = uq_simple_parser(varargin, parse_keys, parse_types);
    % 'beta' option additionally prints beta
    beta_flag = strcmp(uq_cline{1}, 'true');

    % 'theta' option additionally prints beta
    theta_flag = strcmp(uq_cline{2}, 'true');
    
    % 'F' option additionally prints the information matrix
    F_flag = strcmp(uq_cline{3}, 'true');
       
    % 'optim' option additionally prints extensive optimization results
    % report
    optim_flag = strcmp(uq_cline{4}, 'true');
    
    % 'trend' option additionally prints extensive trend-related report
    % report
    trend_flag = strcmp(uq_cline{5}, 'true');
    
    % 'gp' option additionally prints extensive GP-related report
    % report
    gp_flag = any(strcmp({uq_cline{6},uq_cline{7}}, 'true'));
    
    % 'R' option additionally prints R matrix
    R_flag = strcmp(uq_cline{8}, 'true');
    
        
    flagWasSet = beta_flag | theta_flag | F_flag | optim_flag | ...
        trend_flag | gp_flag | R_flag;
else
    flagWasSet = false;
end

%% Produce the fixed header
fprintf('\n%%--------------- Kriging metamodel ---------------%%\n');
fprintf('\tObject Name:\t\t%s\n', KRGModel.Name);

%% If some flag(s) print out the specified elements
if flagWasSet
    
    for ii =  1 : length(outArray)
        current_output = outArray(ii);
        if length(outArray) > 1
            fprintf('--- Output #%i:\n', current_output);
        end
    end
    
    if beta_flag
            fprintf('\n\tRegression Coefficients (beta):\n')
            fprintf('%s\n', ...
                add_leadingChars(uq_sprintf_mat(KRGModel.Kriging(current_output).beta) ,...
                sprintf('\t\t\t\t\t') ));
    end
    
    if F_flag
        fprintf('\nInformation matrix (F):\n')
        fprintf('%s\n', ...
            add_leadingChars(uq_sprintf_mat(KRGModel.Internal.Kriging(current_output).Trend.F) ,...
            sprintf('\t\t\t\t\t') ));
    end
   
    if R_flag
        fprintf('\nCorrelation matrix (R):\n')
        fprintf('%s\n', ...
            add_leadingChars(uq_sprintf_mat(KRGModel.Internal.Kriging(current_output).GP.R) ,...
            sprintf('\t\t\t\t\t') ));
    end
    
    
   if theta_flag
            fprintf('\n\tHyperparameter values (theta):\n')
            fprintf('%s\n', ...
                add_leadingChars(uq_sprintf_mat(KRGModel.Kriging(current_output).theta) ,...
                sprintf('\t\t\t\t\t') ));
    end

    fprintf('%%--------------------------------------------------%%\n');
    return
end


%% Produce the default printout

M = KRGModel.Internal.Runtime.M;
fprintf('\tInput Dimension:\t%i\n', M);
fprintf('\n\tExperimental Design\n')
fprintf('\t\tSampling:\t%s\n', KRGModel.ExpDesign.Sampling)
fprintf('\t\tX size:\t\t[%s]\n', [num2str(size(KRGModel.ExpDesign.X,1)),'x',num2str(size(KRGModel.ExpDesign.X,2))])
fprintf('\t\tY size:\t\t[%s]\n', [num2str(size(KRGModel.ExpDesign.Y,1)),'x',num2str(size(KRGModel.ExpDesign.Y,2))])

for ii =  1 : length(outArray)
    current_output = outArray(ii);
    if length(outArray) > 1
        fprintf('--- Output #%i:\n', current_output);
    end
    %% Trend
    fprintf('\n\tTrend\n')
    fprintf('\t\tType:\t\t%s\n', KRGModel.Internal.Kriging(current_output).Trend.Type)
    if ~any(strcmpi(KRGModel.Internal.Kriging(current_output).Trend.Type, {'simple', 'custom'}))
        fprintf('\t\tDegree:\t\t%i\n', full(KRGModel.Internal.Kriging(current_output).Trend.Degree));
    end
    %% GP
    fprintf('\n\tGaussian Process\n')
    % Check whether a user defined or the default eval_R is used:
    CorrHandle = KRGModel.Internal.Kriging(current_output).GP.Corr.Handle;
    if strcmp(char(CorrHandle), 'uq_eval_Kernel')
        if KRGModel.Internal.Kriging(current_output).GP.Corr.Isotropic
            corrIsotropy = 'isotropic';
        else
            corrIsotropy = 'anisotropic';
        end
        fprintf('\t\tCorr. Type:\t%s(%s)\n', KRGModel.Internal.Kriging(current_output).GP.Corr.Type,corrIsotropy)
        fprintf('\t\tCorr. family:\t%s\n', KRGModel.Internal.Kriging(current_output).GP.Corr.Family)
        fprintf('\t\tsigma^2:\t\t%s\n', KRGModel.Internal.Kriging(current_output).GP.sigmaSQ)
        switch lower(KRGModel.Internal.Kriging(current_output).GP.EstimMethod)
            
            case 'ml'
                fprintf('\tEstimation method:\t%s\n', 'Maximum-Likelihood')
            case 'cv'
                fprintf('\tEstimation method:\t%s\n', 'Cross-Validation')
        end
    else
        % If it is a user-defined handle just show the handle for now
        fprintf('\t\tCorr. Handle:\t%s\n', func2str(CorrHandle))
    end
    %% Optimization
    fprintf('\n\tHyperparameters\n')
    theta = KRGModel.Kriging(current_output).theta ;
    fprintf('\t\ttheta:\t\t[%s]\n', uq_sprintf_mat(theta))
    switch lower(KRGModel.Internal.Kriging(current_output).Optim.Method)
        
        case 'hga'
            fprintf('\t\tOptim. method:\t\t%s\n', 'Hybrid Genetic Algorithm' )
        case 'ga'
            fprintf('\t\tOptim. method:\t\t%s\n', 'Genetic Algorithm' )
        case 'bfgs'
            fprintf('\t\tOptim. method:\t\t%s\n', 'BFGS' )
        case 'sade'
            fprintf('\t\tOptim. method:\t\t%s\n', 'Self-Adaptive Differential Evolution' )
        case 'de'
            fprintf('\t\tOptim. method:\t\t%s\n', 'Differential Evolution' )     
        case 'cmaes'
            fprintf('\t\tOptim. method:\t\t%s\n', 'Covariance matrix adaptation - evolution scheme' ) ;   
       case 'hcmaes'
            fprintf('\t\tOptim. method:\t\t%s\n', 'Hybrid Covariance matrix adaptation - evolution scheme' )
        otherwise
            fprintf('\t\tOptim. method:\t\t%s\n', ...
                KRGModel.Internal.Kriging(current_output).Optim.Method)
    end
    
    
    
    fprintf('\n\tLeave-one-out error:\t%13.7e\n\n', ...
        KRGModel.Error(current_output).LOO);
    if isfield(KRGModel.Error,'Val')
        fprintf('\tValidation error:\t\t%13.7e\n\n', ...
            KRGModel.Error(outArray(ii)).Val);
    end

    
end

fprintf('%%--------------------------------------------------%%\n');


end


function str = add_leadingChars(str, chars , omitLast)

if nargin < 3
    omitLast = 1;
end

str = [chars str] ;

ind = strfind(str, sprintf('\n')) ;

if omitLast
    ind = ind(1:end-1);
end


nCh = length(chars) ;
nCh_tot = 0;
for ii = ind
    str = [str(1 : ii+nCh_tot), chars, str(ii+nCh_tot+1 : end) ] ;
    nCh_tot = nCh_tot + nCh;
end


end



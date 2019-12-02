% name of Matlab file representing the turbine data
turbineName = 'NM80';

%% model description 
% name of Matlab file representing the model
Model.mHandle = @aero_module;
% Optionally, one can pass parameters to the model stores in the cell
% array P
P = getParameterAeroModule(turbineName);
Model.Parameters = P;
Model.isVectorized = false;

%% list of UQ methods to be used for analysis

% specify a list of options from the following list:
methods = {'PCE_LARS'}; % {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};

% for MC, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;
% number of samples with MC
NsamplesMC = [1e1];

% for PCE-Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:4; %[1 2 3 4 5 6];

% % for PCE-OLS:
NsamplesOLS = [8]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat
 
% for PCE-LARS:
NsamplesLARS = [4]; % if not specified, the number of samples from Quad is taken
LARS_repeat = 1; % like MC_repeat

%% Assemble the Input.Marginal for sensitivity analysis by text comparison
ndim = length(P{26});
ntot = length(P{25}.Marginals); 
for i=1:ndim    
    for j = 1:ntot
        if(strcmp([P{25}.Marginals(j).Name,num2str(P{25}.Marginals(j).Index)],[P{26}{i}{1},num2str(P{26}{i}{2})]))
            Input.Marginals(i).Name =  [P{26}{i}{1},num2str(P{26}{i}{2})];
            Input.Marginals(i).Type = P{25}.Marginals(j).Type; 
            Input.Marginals(i).Parameters = P{25}.Marginals(j).Parameters;
            Input.Marginals(i).Bounds = P{25}.Marginals(j).Bounds;
            break;
        end
    end
end
 
% Name of Matlab file representing the turbine data
turbineName = 'NM80'; % 'NM80', 'AVATAR'
% check NM80.m or AVATAR.m or (turbine_name).m for turbine-specific
% settings and definition of uncertainties

%% model description 
% Name of Matlab file representing the model
Model.mHandle = @aero_module;
% Optionally, one can pass parameters to the model stored in the cell
% array P
P = getParameterAeroModule(turbineName);
Model.Parameters = P;
Model.isVectorized = false;

%% list of UQ methods to be used for analysis

% specify a list of options from the following list:

% methods = {'PCE_OLS'}; % {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};

methods = {'PCE_LARS'}; % {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};


% for MC, specify number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;
% number of samples with MC
NsamplesMC = [8 16 32];

% for PCE_Quad, specify the polynomial degrees to be tested
DegreesQuad = 1:3; %[1 2 3 4 5 6];

% % for PCE-OLS:
NsamplesOLS = [8 16 32 64 128]; % if not specified, the number of samples from Quad is taken
OLS_repeat = 1; % like MC_repeat
 
% for PCE-LARS:

NsamplesLARS = [8]; % if not specified, the number of samples from Quad is taken

LARS_repeat = 1; % like MC_repeat

%% Assemble the Input.Marginal for sensitivity analysis by text comparison
% NOTE: check getParameterAeroModule.m to see the definition of the P array
% P{26} contains the uncertain parameters for which we will do sensitivity analysis
ndim = length(P{26}); 
% P{25} contains all possible parameters, deterministic and uncertain, of
% which a subset is used in the sensitivity study (as defined in P{25})
ntot = length(P{25}.Marginals); 
discrete_index = [];
cont_index = [];
discrete_param_vals = [];

% loop over P{26} and for each uncertain parameter get the distribution as
% stored in P{25}
for i=1:ndim    
    for j = 1:ntot
        % find which index we need by looking in struct P{25}
        % store the required information in Input.Marginals(i), which will
        % be used by UQLab
        if(strcmp([P{25}.Marginals(j).Name,num2str(P{25}.Marginals(j).Index)],[P{26}{i}{1},num2str(P{26}{i}{2})]))
            Input.Marginals(i).Name =  [P{26}{i}{1},num2str(P{26}{i}{2})];
            Input.Marginals(i).Type = P{25}.Marginals(j).Type; 
            Input.Marginals(i).Parameters = P{25}.Marginals(j).Parameters;
            Input.Marginals(i).Bounds = P{25}.Marginals(j).Bounds;
            
            if(P{26}{i}{3} ==1) % Get the index and parameter of discrete variable
                discrete_index = [discrete_index i];
                discrete_param_vals = [discrete_param_vals Input.Marginals(i).Parameters(2)];
            else
                cont_index = [cont_index i];
            end
            break;
        end
    end
end



%% model description 

Model.mHandle = @Ishigami;

% name of Matlab file representing the model
% see e.g. Smith, p. 329
a = [1;2];
Model.Parameters = a;

%% input description

% take options from following list:
% methods = {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};
methods = {'MC','PCE_Quad'};

% number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;

DegreesPCE = 1:8; %[1 2 3 4 5 6];
NsamplesMC = [1e1 1e2 1e3 1e4];



%% Description of input marginals
ndim = 3;

%% uniform:
for ii = 1 : ndim
    Input.Marginals(ii).Type = 'Uniform'; 
    Input.Marginals(ii).Parameters = [-pi, pi];
    Input.Marginals(ii).Bounds = [-pi, pi]; % not necessary for uniform but useful for plotting
end

% exact mean and std
mean_exact = a(1)/2;
mean_ref   = mean_exact;
var_exact  = (a(1)^2)/8 + (a(2)*pi^4)/5 + (a(2)^2)*(pi^8)/18 + 1/2;
std_exact  = sqrt(var_exact);
std_ref    = std_exact;

% exact Sobol indices
D1 = (a(2)*pi^4)/5 + (a(2)^2)*(pi^8)/50 + 1/2;
D2 = (a(1)^2)/8;
D3 = 0;
D12 = 0;
D13 = (a(2)^2)*(pi^4)/18 - (a(2)^2)*(pi^8)/50;
D23 = 0;
D123 = 0;

% first order indices
S1 = D1/var_exact;
S2 = D2/var_exact;
S3 = D3/var_exact;
S12 = D12/var_exact;
S13 = D13/var_exact;
S23 = D23/var_exact;
S123 = D123/var_exact;
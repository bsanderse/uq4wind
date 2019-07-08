%% model description 

Model.mHandle = @linear_portfolio;

% name of Matlab file representing the model
% see e.g. Smith, p. 321, 328
a = [2;1];
Model.Parameters = a;

%% input description

% take options from following list:
% methods = {'MC','PCE_Quad','PCE_OLS','PCE_LARS'};
methods = {'MC','PCE_Quad'};

% number of times to repeat MC-based methods to obtain 'nice' convergence
% graphs
MC_repeat = 1;

DegreesPCE = 1:6; %[1 2 3 4 5 6];
NsamplesMC = [1e1 1e2 1e3 1e4];



%% Description of input marginals
ndim = 2;

%% uniform:
% for ii = 1 : ndim
%     Input.Marginals(ii).Type = 'Uniform'; 
%     Input.Marginals(ii).Parameters = [0, 1];
%     Input.Marginals(ii).Bounds = [0, 1]; % not necessary for uniform but useful for plotting
% end
% 
% % exact mean:
% mean_exact = 0;
% 
% for i=1:ndim
%     
%     mean_exact = mean_exact + 0.5*a(i)*(Input.Marginals(i).Bounds(2)^2 - Input.Marginals(i).Bounds(1)^2);
%     
% end

%% Gaussian 
mu    = [0;0];
sigma = [1;3];
for ii = 1 : ndim
    Input.Marginals(ii).Type = 'Gaussian'; 
    Input.Marginals(ii).Parameters = [mu(ii), sigma(ii)];
    % don't use Bounds here, as it will truncate the Gaussian
end

% exact mean:
% mean_exact = integral2( @(x,y) 1/(2*pi) * (a(1)*x+a(2)*y).*exp(-x.^2/2-y.^2/2),-Inf,Inf,-Inf,Inf);
% var_exact  = integral2( @(x,y) 1/(2*pi) * ((a(1)*x+a(2)*y).^2).*exp(-x.^2/2-y.^2/2),-Inf,Inf,-Inf,Inf) - mean_exact^2
mean_exact = 0;
var_exact  = a(1)^2*sigma(1)^2 + a(2)^2*sigma(2)^2;
std_exact  = sqrt(var_exact);
mean_ref   = 1; % because mean_exact = 0
std_ref    = std_exact;


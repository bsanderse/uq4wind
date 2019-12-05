%% simple UQLab test with Ishigami function
% see also https://www.uqlab.com/model-parameters

clearvars
close all

%% start uqlab
uqlab

%% model and input specification

% specify file with matlab model
Model.mFile = 'Ishigami'; % or: Model.mHandle = @Ishigami
% model parameters
a     = [7;0.1];
Model.Parameters = a;

% create and add the model to UQLab
myModel = uq_createModel(Model);

% create 3 input parameters ('X') and specify distribution
for ii = 1 : 3
    Input.Marginals(ii).Type = 'Uniform'; 
    Input.Marginals(ii).Parameters = [-pi, pi];
end

% create input object in UQLab
myInput = uq_createInput(Input);

% display input properties
uq_print(myInput);
uq_display(myInput);


%% let's now perform some UQ analysis

%% Monte Carlo: 
disp('==============================================================')
disp('Monte Carlo sampling')
% number of samples
Nsamples = 100;
% get random sample of the inputs ('experimental design')
X_ED = uq_getSample(Nsamples, 'MC');
% evaluate the model at these inputs
Y_ED = uq_evalModel(myModel,X_ED);
% get mean and standard deviation            
mean_MC = mean(Y_ED)
std_MC  = std(Y_ED)

% make plots of the response:
% uq_figure
% uq_histogram(Y_ED)


%% PCE - OLS:
disp('==============================================================')
disp('Polynomial Chaos Expansion')
% number of samples
NsamplesLARS = 100;

% UQLab settings to set up PCE with ordinary least squares (OLS)
metamodelLARS.FullModel = myModel;
metamodelLARS.Input     = myInput;
metamodelLARS.Type      = 'Metamodel';
metamodelLARS.MetaType  = 'PCE';
metamodelLARS.Method    = 'LARS';

% specify array of possible degrees;
% the degree with the lowest Leave-One-Out cross-validation error (LOO error)
% is automatically selected:
metamodelLARS.Degree    = 2:10;
metamodelLARS.TruncOptions.qNorm = 0.75;

metamodelLARS.ExpDesign.NSamples = NsamplesLARS;
metamodelLARS.ExpDesign.Sampling = 'LHS'; % latin hypercube sampling
myPCE_LARS = uq_createModel(metamodelLARS);

% moments of solution            
mean_LARS = myPCE_LARS.PCE.Moments.Mean
std_LARS  = sqrt(myPCE_LARS.PCE.Moments.Var)     


%%
disp('==============================================================')
disp('Exact results')
% exact mean and std
mean_exact = a(1)/2
var_exact  = (a(1)^2)/8 + (a(2)*pi^4)/5 + (a(2)^2)*(pi^8)/18 + 1/2;
std_exact  = sqrt(var_exact)
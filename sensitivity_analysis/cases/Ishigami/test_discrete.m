clear all
clc

% add path
% UQLab_path = '/Users/sanderse/Dropbox/work/Programming/UQ/UQLabCore_Rel1.0.0/';
% run('../../config.m');
% addpath(genpath(UQLab_path));

% start uqlab
uqlab

% specify file with matlab model
Model.mFile = 'Ishigami'; % or: Model.mHandle = @Ishigami
% model parameters
a     = [7;0.1];
Model.Parameters = a;

% create and add the model to UQLab
myModel = uq_createModel(Model);



% uniform:
ncont=2;
ndisc=1;

for ii = 1 : ncont
    Input.Marginals(ii).Type = 'Uniform';
    Input.Marginals(ii).Parameters = [-pi, pi];
    Input.Marginals(ii).Bounds = [-pi, pi]; % not necessary for uniform but useful for plotting
end
Input.Marginals(3).Type = 'Constant';
Input.Marginals(3).Parameters = 2;
    
    
NsamplesLARS = 100;

% create input object with UQLab
myInput = uq_createInput(Input);


metamodelLARS.FullModel = myModel;
metamodelLARS.Input     = myInput;
metamodelLARS.Type      = 'Metamodel';
metamodelLARS.MetaType  = 'PCE';
metamodelLARS.Method    = 'LARS';
metamodelLARS.Degree    = 1:4; % this automatically switches on degree adaptive PCE
metamodelLARS.TruncOptions.qNorm = 0.75;


% use manual experimental design:
%         X_ED = uq_getSample(NsamplesLARS(i),'MC') ;
%         Y_ED = uq_evalModel(myModel,X_ED);
%         metamodelLARS.ExpDesign.X = X_ED;
%         metamodelLARS.ExpDesign.Y = Y_ED;


% use sampling strategy, note that default is MC!
metamodelLARS.ExpDesign.Sampling = 'LHS'; % or 'LHS' or 'Sobol' or 'Halton'
metamodelLARS.ExpDesign.NSamples = NsamplesLARS;

X_ED = uq_getSample(myInput,NsamplesLARS,'LHS');

Xnew = repmat(X_ED,2,1);
Xnew(101:end,3) = 3;

Ynew = uq_evalModel(myModel,Xnew);

metamodelLARS.ExpDesign.X = Xnew;
metamodelLARS.ExpDesign.Y = Ynew;

metamodelLARS.ExpDesign.Sampling = 'user'; % or 'LHS' or 'Sobol' or 'Halton'
metamodelLARS.ExpDesign.NSamples = 200;

myPCE_LARS     = uq_createModel(metamodelLARS);

% X_ED = myPCE_LARS.ExpDesign.X;

% moments of solution
% mean_LARS(k,i) = myPCE_LARS.PCE.Moments.Mean;
% std_LARS(k,i)  = sqrt(myPCE_LARS.PCE.Moments.Var);
% 
% % Sobol analysis
% % note the same options structure SobolOpts can be re-used to create a new analysis on the PCE model
% SobolAnalysis_LARS    = uq_createAnalysis(SobolOpts);
% SobolResults_LARS     = SobolAnalysis_LARS.Results;
% Sobol_LARS_FirstOrder(k, i, 1:ncont) = SobolResults_LARS.FirstOrder;
% Sobol_LARS_Total(k, i, 1:ncont)      = SobolResults_LARS.Total;

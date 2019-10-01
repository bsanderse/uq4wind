function pass = uq_PCE_test_qNormAdapt( level )
% PASS = UQ_PCE_TEST_QNORMADAPT(LEVEL): running tests for qNorm
% adaptive PCE and check if they get the 'correct' qNorm.

% Initialize test:
pass = 1;
evalc('uqlab');
rng(500);
if nargin < 1
    level = 'normal';
end
fprintf(['\nRunning: |' level '| uq_PCE_test_adaptqNormMaxInter...\n']);

%% INPUT
% values taken from the default phimecasoft example
% Define the probabilistic model.
for i = 1:3
    Input.Marginals(i).Type = 'Uniform';
    Input.Marginals(i).Parameters = [-pi, pi];
end

myInput = uq_createInput(Input,'-private');

%% MODEL
% Physical model: Ishigami function
modelopts.Name = 'ishigami test';
modelopts.mFile = 'uq_ishigami';
FullModel = uq_createModel(modelopts,'-private');

%% PCE Metamodel
samplesizes = 10:20:200;
qnorms = [    0.5000    0.5000    0.6000    0.9500    0.7000    0.7500    0.7000    0.7000    0.7000    0.6500];

for ii = 1 : length(samplesizes)
    rng(500);

    clear metaopts
    
    metaopts.Type = 'metamodel';
    metaopts.MetaType = 'PCE';
    metaopts.Method = 'LARS';
    metaopts.Input = myInput;
    metaopts.FullModel = FullModel;
    metaopts.qNormEarlyStop = true;
    metaopts.DegreeEarlyStop = true;
    metaopts.LARS.LarsEarlyStop = true;
    
    metaopts.ExpDesign.NSamples = samplesizes(ii);
    
    metaopts.Degree = 5:10;
    metaopts.TruncOptions.qNorm = 0.5:0.05:1;
    
    myPCE = uq_createModel(metaopts,'-private');
    
    %% Validation
    % check the multi indices
    NofVars = sum(myPCE.PCE.Basis.Indices~=0,2);
    pass = pass & myPCE.Internal.PCE.Basis.Truncation.qNorm == qnorms(ii);
end

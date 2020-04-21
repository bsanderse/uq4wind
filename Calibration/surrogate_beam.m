clearvars
rng(100,'twister')
uqlab

%%
uq_figure
[I,~] = imread('SimplySupportedBeam.png');
image(I)
axis equal
set(gca, 'visible', 'off')

%%
ModelOpts.mFile = 'uq_SimplySupportedBeam';
ModelOpts.isVectorized = true;

myForwardModel = uq_createModel(ModelOpts);

%%
PriorOpts.Marginals(1).Name = 'b';               % beam width
PriorOpts.Marginals(1).Type = 'Constant';
PriorOpts.Marginals(1).Parameters = [0.15];      % (m)

PriorOpts.Marginals(2).Name = 'h';               % beam height
PriorOpts.Marginals(2).Type = 'Constant';
PriorOpts.Marginals(2).Parameters = [0.3];       % (m)

PriorOpts.Marginals(3).Name = 'L';               % beam length
PriorOpts.Marginals(3).Type = 'Constant';
PriorOpts.Marginals(3).Parameters = 5;           % (m)

PriorOpts.Marginals(4).Name = 'E';               % Young's modulus
PriorOpts.Marginals(4).Type = 'LogNormal';
PriorOpts.Marginals(4).Moments = [30000 4500];   % (MPa)

PriorOpts.Marginals(5).Name = 'p';               % uniform load
PriorOpts.Marginals(5).Type = 'Gaussian';
PriorOpts.Marginals(5).Moments = [0.012 0.012*0.05]; % (kN/m)

myPriorDist = uq_createInput(PriorOpts);

%%
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.ExpDesign.NSamples = 50;
mySurrogateModel = uq_createModel(MetaOpts);

%%
myData.y = [12.84; 13.12; 12.13; 12.19; 12.67]/1000;  % (m)
myData.Name = 'Mid-span deflection';

%%
BayesOpts.Type = 'Inversion';
BayesOpts.Data = myData;

%%
BayesOpts.ForwardModel.Model = mySurrogateModel;

%%
myBayesianAnalysis_surrogateModel = uq_createAnalysis(BayesOpts);

%%
uq_print(myBayesianAnalysis_surrogateModel)
uq_display(myBayesianAnalysis_surrogateModel)

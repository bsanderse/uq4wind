%% SENSITIVITY ANALYSIS: DEPENDENT INPUT VARIABLES
%
% In this example, the Kucherenko and ANCOVA sensitivity indices are
% computed and compared for the
% <https://www.sfu.ca/~ssurjano/shortcol short column function> with
% correlated input parameters.

%% 1 - INITIALIZE UQLAB
%
% Clear all variables from the workspace, set the random number generator
% for reproducible results, and initialize the UQLab framework:
clearvars
rng(101,'twister')
uqlab

%% 2 - COMPUTATIONAL MODEL
%
% The short-column function is an analytical 3-dimensional function
% that models the limit state of a short steel column with a rectangular
% cross-section, subjected to a bending moment and an axial force.
% The short column function |uq_shortcol| is supplied with UQLab.
%
% Create a MODEL object from the function file:
ModelOpts.mFile = 'uq_shortcol';
myModel = uq_createModel(ModelOpts);

%%
% Type |help uq_shortcol| for information on the model structure.

%% 3 - PROBABILISTIC INPUT MODEL
%
% The probabilistic input model consists of three random variables with the
% following marginals:
%%
% <html>
% <table border=1><tr>
% <td><b>Variable</b></td>
% <td><b>Description</b></td>
% <td><b>Distribution</b></td>
% <td><b>Mean</b></td>
% <td><b>Std. deviation</b></td></tr>
% <tr>
% <td>Y</td>
% <td>Yield stress (MPa)</td>
% <td>Lognormal</td>
% <td>5</td>
% <td>0.5</td>
% </tr>
% <tr>
% <td>M</td>
% <td>Bending moment (N.mm)</td>
% <td>Gaussian</td>
% <td>2000</td>
% <td>400</td>
% </tr>
% <tr>
% <td>P</td>
% <td>Axial force (N)</td>
% <td>Gaussian</td>
% <td>500</td>
% <td>100</td>
% </tr>
% </table>
% </html>

%%
% Define the marginal distributions of the input model:
InputOpts.Marginals(1).Name = 'Y'; % yield stress
InputOpts.Marginals(1).Type = 'Lognormal';
InputOpts.Marginals(1).Moments = [5 0.5]; % (MPa)

InputOpts.Marginals(2).Name = 'M'; % bending moment
InputOpts.Marginals(2).Type = 'Gaussian';
InputOpts.Marginals(2).Moments = [2000 400]; % (N.mm)

InputOpts.Marginals(3).Name = 'P'; % axial force
InputOpts.Marginals(3).Type = 'Gaussian';
InputOpts.Marginals(3).Moments = [500 100]; % (N)

%%
% The parameters $P$ and $M$ are assumed to be correlated
% under Gaussian copula with a rank correlation of 0.72.
%
% Define the Gaussian copula and its parameters:
InputOpts.Copula.Type = 'Gaussian';
InputOpts.Copula.RankCorr = ...
    [1        0         0      ;...
    0         1         0.72;...
    0         0.72   1      ];

%%
% Create an INPUT object based on the marginals and copula:
myInput = uq_createInput(InputOpts);

%%
% Show the correlated input:
uq_display(myInput)

%% 4 - SENSITIVITY ANALYSIS
%
% In the following, the Kucherenko analysis is compared
% to the ANCOVA indices.

%% 4.1 Kucherenko indices
% 
% First, select the sensitivity tool in UQLab
% and specify Kucherenko analysis:
KucherenkoOpts.Type = 'Sensitivity';
KucherenkoOpts.Method = 'Kucherenko';

%%
% Specify the sample size for the analysis (the resulting cost depends 
% on the estimator):
KucherenkoOpts.Kucherenko.SampleSize = 5e4;

%%
% Run the Kucherenko analysis:
KucherenkoAnalysis = uq_createAnalysis(KucherenkoOpts);

%% 4.2 ANCOVA indices
%
% Select the sensitivity tool in UQLab and specify ANCOVA analysis:
ANCOVAOpts.Type = 'Sensitivity';
ANCOVAOpts.Method = 'ANCOVA';

%%
% Specify the sample size for the analysis
% (it will be used to construct the PCE):
ANCOVAOpts.ANCOVA.SampleSize = 150;

%%
% Run the ANCOVA sensitivity analysis:
ANCOVAAnalysis = uq_createAnalysis(ANCOVAOpts);

%% 5 - RESULTS
%
% Print out a report on the Kucherenko analysis:
uq_print(KucherenkoAnalysis)

%%
% Display a graphical representation of the analysis:
uq_display(KucherenkoAnalysis)

%%
% Print out a report on the ANCOVA analysis:
uq_print(ANCOVAAnalysis)

%%
% Display a graphical representation of the analysis:
uq_display(ANCOVAAnalysis)
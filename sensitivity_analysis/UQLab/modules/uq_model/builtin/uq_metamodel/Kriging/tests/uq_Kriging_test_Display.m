function success = uq_Kriging_test_Display(level)
%UQ_KRIGING_TEST_DISPLAY tests the uq_display on Kriging objects.
%
%   Summary:
%   Make sure that uq_display works on a Kriging metamodel objects. The
%   test if uq_display throws error, but it does not check for correctness
%   in appearance of the resulting figures.
%
%   SUCCESS = UQ_KRIGING_TEST_DISPLAY(LEVEL) carried out non-regression
%   tests with the test depth specified in the string LEVEL for the UQLab
%   display feature (via uq_display function) of Kriging metamodel objects.

%% Initialize UQLab
uqlab('-nosplash')
rng(423,'twister')

if nargin < 1
    level = 'normal';
end
fprintf('\nRunning: | %s | uq_Kriging_test_Display...\n',level)

success = true;

%% Define common options

% Kriging metamodel
MetaOpts.Type =  'Metamodel';
MetaOpts.MetaType = 'Kriging';
MetaOpts.Display = 'quiet';

%% Problem 1 Setup: One-Dimensional Function, Scalar Output (Interpolation)
X = -pi + 2*pi*rand(10,1);
Y = X .* sin(X);

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Create the Kriging metamodel
myKriging = uq_createModel(MetaOpts,'-private');

% Create the plots
try
    uq_display(myKriging)
    close(gcf)
    uq_display(myKriging, 1, 'R')
    close(gcf)
catch e
    rethrow(e)
end

%% Problem 2 Setup: Two-Dimensional Function, Scalar Output (Interpolation)

% Branin-Hoo function
% Create a UQLab INPUT object:
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [-5 10];

InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0 15];

myInput = uq_createInput(InputOpts,'-private');

X = uq_getSample(myInput,50);
Y = uq_branin(X);

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Create the Kriging metamodel
myKriging = uq_createModel(MetaOpts,'-private');

% Create the plots
try
    uq_display(myKriging)
    close(gcf)
    close(gcf)
    uq_display(myKriging, 1, 'R')
    close(gcf)
catch e
    rethrow(e)
end

%% Problem 3 Setup: One-Dimensional Function, Vector Output (Interpolation)
X = -pi + 2*pi*rand(10,1);
Y1 = X .* sin(X);
Y2 = 1 + X*5e-2 + sin(X) ./ X;
Y = [Y1 Y2];

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Create the Kriging metamodel
myKriging = uq_createModel(MetaOpts,'-private');

% Create the plots
try
    uq_display(myKriging)
    close(gcf)
    uq_display(myKriging,1:2)
    close(gcf)
    close(gcf)
    uq_display(myKriging,[2 1])
    close(gcf)
    close(gcf)
    uq_display(myKriging,1)
    close(gcf)
    uq_display(myKriging,2)
    close(gcf)
    uq_display(myKriging, 1, 'R')
    close(gcf)
    uq_display(myKriging, 2, 'R')
    close(gcf)
    uq_display(myKriging, 1:2, 'R')
    close(gcf)
    close(gcf)
    uq_display(myKriging, [2 1], 'R')
    close(gcf)
    close(gcf)
catch e
    rethrow(e)
end

%% Problem 4 Setup: Two-Dimensional Function, Vector Output (Interpolation)

% Branin-Hoo function
% Create a UQLab INPUT object:
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [-5 10];

InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0 15];

myInput = uq_createInput(InputOpts,'-private');

X = uq_getSample(myInput,50);
Y = [uq_branin(X) uq_branin(X)];

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Create the Kriging metamodel
myKriging = uq_createModel(MetaOpts,'-private');

% Create the plots
try
    uq_display(myKriging)
    close(gcf)
    close(gcf)
    uq_display(myKriging,1:2)
    close(gcf)
    close(gcf)
    close(gcf)
    close(gcf)
    uq_display(myKriging,[2 1])
    close(gcf)
    close(gcf)
    close(gcf)
    close(gcf)
    uq_display(myKriging,1)
    close(gcf)
    close(gcf)
    uq_display(myKriging,2)
    close(gcf)
    close(gcf)
    uq_display(myKriging, 1, 'R')
    close(gcf)
    uq_display(myKriging, 2, 'R')
    close(gcf)
    uq_display(myKriging, 1:2, 'R')
    close(gcf)
    close(gcf)
    uq_display(myKriging, [2 1], 'R')
    close(gcf)
    close(gcf)
catch e
    rethrow(e)
end

%% Problem 5 Setup: One-Dimensional Function, Scalar Output (Regression)
X = -pi + 2*pi*rand(50,1);
Y = X .* sin(X) + 1.0*randn(size(X,1),1);

% Regression Options
SigmaNSQs = {'auto', true, 1.0, eye(size(X,1))};

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

for SigmaNSQ = SigmaNSQs
    try
        MetaOpts.Regression.SigmaNSQ = SigmaNSQ{1};
        % Create the Kriging metamodel
        myKriging = uq_createModel(MetaOpts,'-private');
        % Create the plots
        uq_display(myKriging)
        close(gcf)
        uq_display(myKriging, 1, 'R')
        close(gcf)
    catch e
        rethrow(e)
    end
end

%% Problem 6 Setup: Two-Dimensional Function, Scalar Output (Regression)

% Branin-Hoo function
% Create a UQLab INPUT object:
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [-5 10];

InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0 15];

myInput = uq_createInput(InputOpts,'-private');

X = uq_getSample(myInput,50);
Y = uq_branin(X) + 1.0*randn(50,1);

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Create the Kriging metamodel
myKriging = uq_createModel(MetaOpts,'-private');

% Create the plots
try
    uq_display(myKriging)
    close(gcf)
    close(gcf)
    uq_display(myKriging, 1, 'R')
    close(gcf)
catch e
    rethrow(e)
end

%% Problem 7 Setup: One-Dimensional Function, Vector Output (Regression)
X = -pi + 2*pi*rand(50,1);
Y = X .* sin(X) + 1.0*randn(size(X,1),1);
Y = [Y Y];

% Regression Options
SigmaNSQs = {'auto', true, 1.0, eye(size(X,1))};

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

for SigmaNSQ = SigmaNSQs
    try
        MetaOpts.Regression.SigmaNSQ = SigmaNSQ{1};
        % Create the Kriging metamodel
        myKriging = uq_createModel(MetaOpts,'-private');
        % Create the plots
        uq_display(myKriging)
        close(gcf)
        uq_display(myKriging,1:2)
        close(gcf)
        close(gcf)
        close(gcf)
        close(gcf)
        uq_display(myKriging,[2 1])
        close(gcf)
        close(gcf)
        close(gcf)
        close(gcf)
        uq_display(myKriging,1)
        close(gcf)
        uq_display(myKriging,2)
        close(gcf)
        uq_display(myKriging, 1, 'R')
        close(gcf)
        uq_display(myKriging, 2, 'R')
        close(gcf)
        uq_display(myKriging, 1:2, 'R')
        close(gcf)
        close(gcf)
        uq_display(myKriging, [2 1], 'R')
        close(gcf)
        close(gcf)
    catch e
        rethrow(e)
    end
end

%% Problem 8 Setup: Two-Dimensional Function, Vector Output (Regression)

% Branin-Hoo function
% Create a UQLab INPUT object:
InputOpts.Marginals(1).Type = 'Uniform';
InputOpts.Marginals(1).Parameters = [-5 10];

InputOpts.Marginals(2).Type = 'Uniform';
InputOpts.Marginals(2).Parameters = [0 15];

myInput = uq_createInput(InputOpts,'-private');

X = uq_getSample(myInput,50);
Y = [uq_branin(X)+randn(50,1) uq_branin(X)+randn(50,1)];

% Experimental design
MetaOpts.ExpDesign.X = X;
MetaOpts.ExpDesign.Y = Y;

% Create the Kriging metamodel
myKriging = uq_createModel(MetaOpts,'-private');

% Create the plots
try
    uq_display(myKriging)
    close(gcf)
    close(gcf)
    uq_display(myKriging,1:2)
    close(gcf)
    close(gcf)
    close(gcf)
    close(gcf)
    uq_display(myKriging,[2 1])
    close(gcf)
    close(gcf)
    close(gcf)
    close(gcf)
    uq_display(myKriging,1)
    close(gcf)
    close(gcf)
    uq_display(myKriging,2)
    close(gcf)
    close(gcf)
    uq_display(myKriging, 1, 'R')
    close(gcf)
    uq_display(myKriging, 2, 'R')
    close(gcf)
    uq_display(myKriging, 1:2, 'R')
    close(gcf)
    close(gcf)
    uq_display(myKriging, [2 1], 'R')
    close(gcf)
    close(gcf)
catch e
    rethrow(e)
end

end
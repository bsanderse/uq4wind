function pass = uq_Kriging_selftest(level)
% PASS = UQ_KRIGING_SELFTEST(LEVEL): suite of non-regression and consistency
%     checks for the Kriging module of UQLab
%
% See also: uq_selftest_uq_metamodel

%% initialize test
uqlab('-nosplash');

if nargin < 1
    level = 'normal'; % TBD: Time that the tests will take
end   

pass = 0;

%% Are the required toolboxes available?
% Required toolboxes for the Kriging module
req_toolbox_names = {'Optimization Toolbox', 'Global Optimization Toolbox'};
% Check 
[ret_checks, ret_names] = uq_check_toolboxes();
OPTIM_TOOLBOX_EXISTS = any(strcmpi(req_toolbox_names{1}, ret_names(ret_checks)));
GOPTIM_TOOLBOX_EXISTS = any(strcmpi(req_toolbox_names{2}, ret_names(ret_checks)));

%% Test Names are defined here
TestNames = {'uq_Kriging_test_TrendTypes'};

if OPTIM_TOOLBOX_EXISTS 
    TestNames = [TestNames, {'uq_Kriging_test_ExpDesigns'}];
else
    warning('Skipping some tests due to absense of the Optimization Toolbox')
end

if OPTIM_TOOLBOX_EXISTS && GOPTIM_TOOLBOX_EXISTS
    TestNames = [TestNames, {'uq_Kriging_test_OptimResult', ...
        'uq_Kriging_test_CustomKriging','uq_Kriging_test_Constant'}];
else
    warning('Skipping some tests due to absense of the Global Optimization Toolbox and/or the Optimization Toolbox')
end
          
%% Recursively execute each test defined in TestNames
success = zeros(length(TestNames),1);
Times = zeros(length(TestNames),1);
TestTimer = tic;
Tprev = 0;
for iTest = 1 : length(TestNames)
    % obtain the function handle of current test from its name
    testFuncHandle = str2func(TestNames{iTest});
    % run test
    success(iTest) = testFuncHandle(level);
    % calculate the time required from the current test to execute
    Times(iTest) = toc(TestTimer) - Tprev ;
    Tprev = Times(iTest);
end


%% Print out the results table and info:
Result = {'ERROR','OK'};
ResultChar = 60; % Character where the result of test is displayed
MinusLine(1:ResultChar+7) = deal('-');
fprintf('\n%s\n',MinusLine);
fprintf('UQ_SELFTEST_UQ_KRIGING RESULTS');
fprintf('\n%s\n',MinusLine);
for ii = 1:length(success)
    points(1:max(2,ResultChar-size(TestNames{ii},2))) = deal('.');
    fprintf('%s %s %s @ %g sec.\n',TestNames{ii},points,Result{success(ii)+1},Times(ii));
    clear points
end
fprintf('%s\n',MinusLine);

%% Where all tests passed?  If not final pass = 0
if all(success)
    pass = 1;
    fprintf('\n');
    fprintf(['SUCCESS: uq_Kriging module ' level ' test was successful.\n']);
else
    
end
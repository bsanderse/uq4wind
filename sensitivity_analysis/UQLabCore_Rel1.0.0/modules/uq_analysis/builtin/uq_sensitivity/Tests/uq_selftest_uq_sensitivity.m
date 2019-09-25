function pass = uq_selftest_uq_sensitivity(level)
% PASS = UQ_SELFTEST_UQ_SENSITIVITY(LEVEL): non-regression and consistency
%     test-suite for the Sensitivity module of UQLab. A number of tests are
%     executed and their results are summarized and printed on screen.
%
% See also: UQ_SENSITIVITY

uqlab('-nosplash');
if nargin < 1
    level = 'normal';
end

% Tests to perform:
TestNames = {'uq_test_correlation',...
    'uq_test_correlation_constant',...
    'uq_test_SRC',...
    'uq_test_perturbation',...
    'uq_test_cotter_basic',...
    'uq_test_sobol_indices',...
    'uq_test_sobol_indices_constant',...
    'uq_test_sobol_indices_PCE',...
    'uq_test_sobol_indices_PCE_constant',...
    'uq_test_morris',...
    'uq_test_morris_high_order_interactions', ...
    'uq_test_morris_linear_model', ...
    'uq_test_morris_trajectory', ...
    'uq_test_sobol_high_order_interactions', ...
    'uq_test_cotter_high_order_interactions', ...
    'uq_test_sensitivity_ishigami_outputs',...
    'uq_test_sobol_indices_LRA',...
    'uq_test_borgonovo_indices',...
    'uq_test_ancova_indices',...
    'uq_test_kucherenko_indices'...
    };

Ntests = length(TestNames);
success = false(1, Ntests);
Times = zeros(1, Ntests);
for ii = 1:Ntests
    TestTimer = tic;
    success(ii) = eval([TestNames{ii} '(level);']);
    Times(ii) = toc(TestTimer);
end

% Print out the results table and info:
Result = {'ERROR','OK'};
ResultChar = 60; % Character where the result of test is displayed
MinusLine(1:ResultChar + 7) = deal('-');
fprintf('\n%s\n',MinusLine);
fprintf('UQ_SELFTEST_UQ_SENSITIVITY RESULTS');
fprintf('\n%s\n',MinusLine);
for ii = 1:Ntests
    points(1:max(2,ResultChar-size(TestNames{ii},2))) = deal('.');
    fprintf('%s %s %s @ %g sec.\n',TestNames{ii},points,Result{success(ii)+1},Times(ii));
    clear points
end
fprintf('%s\n',MinusLine);

if all(success)
    pass = 1;
    fprintf('\n');
    fprintf(['SUCCESS: uq_sensitivity module ' level ' test was successful.\n']);
else
    pass = 0;
    fprintf('\n');
    fprintf(['FAIL: uq_sensitivity module ' level ' test failed.\n']);
end
fprintf('Total time: %g',sum(Times));
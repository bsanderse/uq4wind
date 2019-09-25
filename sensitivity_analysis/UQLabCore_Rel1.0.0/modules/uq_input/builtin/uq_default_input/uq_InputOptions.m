% UQ_INPUTOPTIONS displays a helper for the main options needed to create an 
% probabilistic input in UQLab.
%    UQ_INPUTOPTIONS displays the main options needed by the command
%    <a href="matlab:help uq_createInput">uq_createInput</a> to create an INPUT object in UQLab 
%
%    See also: uq_createInput, uq_getSample, uq_getInput, uq_listInputs,
%              uq_selectInput 
disp('Quickstart guide to the UQLab Input Module')
disp('  ')
disp('In the UQLab software, INPUT objects are created by the command:')
disp('    myInput = uq_createInput(INPUTOPTIONS)')
disp('The options are specified in the INPUTOPTIONS structure.')
disp(' ');
disp('Example: to create a 2-dimensional INPUT object that represents two random variables ')
disp('in [-1,1] with covariance matrix C = [1 0.5; 0.5 1], type:')
disp('    INPUTOPTIONS.Marginals(1).Type = ''Uniform'';')
disp('    INPUTOPTIONS.Marginals(1).Parameters = [-1,1];')
disp('    INPUTOPTIONS.Marginals(2).Type = ''Uniform'';')
disp('    INPUTOPTIONS.Marginals(2).Parameters = [-1,1];')
disp('    INPUTOPTIONS.Copula.Type = ''Gaussian'';')
disp('    INPUTOPTIONS.Copula.Parameters = [1 0.5; 0.5 1];')
disp('    myInput = uq_createInput(INPUTOPTIONS);')
disp(' ');
disp('To draw 100 samples from the resulting object, type:');
disp('    X = uq_getSample(100)');
disp(' ');
AvailableMarginals = uq_getAvailableMarginals;
disp('The following marginal distributions are available for use in UQLab:')
for ii = 1:length(AvailableMarginals)
fprintf('    - %s\n',regexprep(AvailableMarginals{ii},'(\<[a-z])','${upper($1)}'));
end
disp(' ');
disp('Please refer to the Input User Manual (<a href="matlab:uq_doc(''Input'',''html'')">HTML</a>,<a href="matlab:uq_doc(''Input'',''pdf'')">PDF</a>)');
disp('for more detailed information on the available features.')
disp(' ');




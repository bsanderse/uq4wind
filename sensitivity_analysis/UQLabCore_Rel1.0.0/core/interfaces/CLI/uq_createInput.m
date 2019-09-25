% UQ_CREATEINPUT   create a UQLab INPUT object 
%    myInput = UQ_CREATEINPUT(INPUTOPTS) creates and stores in the current
%    UQLab session an object of type INPUT based on the configuration
%    options given in the INPUTOPTS structure. 
%
%    myInput = UQ_CREATEINPUT(INPUTOPTS, '-private') creates an object of
%    type INPUT based on the configuration options given in the INPUTOPTS
%    structure, without saving it in the UQLab session. This is useful,
%    e.g., in case of parametric studies which create a large number of
%    temporary INPUT objects that do not need to be stored in memory.
%
%    Example: create a 2-dimensional INPUT object with independent
%    uniform random variables in [-1,1]:
%       INPUTOPTS.Marginals(1).Type = 'Uniform';
%       INPUTOPTS.Marginals(1).Parameters = [-1 1];
%       INPUTOPTS.Marginals(2).Type = 'Uniform';
%       INPUTOPTS.Marginals(2).Parameters = [-1 1];
%       myInput = UQ_CREATEINPUT(INPUTOPTS);
%
%    To draw N samples from the resulting object use:
%       X = uq_getSample(N)
%
%    For additional information about the available INPUTOPTS configuration
%    fields type: 
%       help <a href="matlab:uq_InputOptions">uq_InputOptions</a>
%
%    To open the INPUT User Manual, type:
%       uq_doc('Input')
%
%    See also: uq_getInput, uq_listInputs, uq_selectInput, uq_getSample,
%              uq_createModel, uq_createAnalysis 
%

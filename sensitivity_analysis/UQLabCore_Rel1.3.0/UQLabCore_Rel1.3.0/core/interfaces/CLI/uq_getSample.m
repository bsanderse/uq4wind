% UQ_GETSAMPLE get a sample from the current INPUT object
%    X = UQ_GETSAMPLE(N) draws a sample of size N from the currently
%    selected UQLab INPUT object. Note that size(X) = [N M], where M is the
%    dimension of the input space.
%
%    X = UQ_GETSAMPLE(myInput,N) draws a sample of size N from the INPUT object
%    myInput  
%
%    X = UQ_GETSAMPLE(..., METHOD) use the specified METHOD to draw samples
%    from the input objects. 
%    METHOD is one of the following strings:
%    'MC'     - Standard Monte Carlo sampling
%    'LHS'    - Latin hypercube sampling (space filling)
%    'Sobol'  - Sobol' pseudorandom sequence
%    'Halton' - Halton pseudorandom sequence
%
%    See also: uq_createInput, uq_selectInput, uq_getInput
%
%    For a list of all available sampling methdos and options, please
%    consult the INPUT User Manual (<a href="matlab:uq_doc('input','html')">HTML</a>,<a href="matlab:uq_doc('input','pdf')">PDF</a>)

% UQ_EVALMODEL evaluate a UQLab MODEL object.
%    Y = UQ_EVALMODEL(X) evaluates the currently selected UQLab MODEL on the
%    vector of input parameters X. Note that size(X) = [N M], where N is
%    the number of realizations of the input parameters, and M is the
%    dimension of the input parameter space.  
%
%    Y = UQ_EVALMODEL(myModel,X) evaluates the UQLab MODEL object myModel
%    on the vector of input parameters X. 
%
%    [Y1,...,YM] = UQ_EVALMODEL(...) returns multiple outputs if the
%    MODEL object supports multiple return values (e.g. Kriging metamodels)
%
%    For more detailed usage scenarios of UQ_EVALMODEL, please refer to the
%    relevant <a href="matlab:uq_doc">UQLab User Manuals</a>
%
%    See also: uq_createModel, uq_getModel, uq_listModels, uq_selectModel
%

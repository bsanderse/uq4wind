% UQ_GETMODEL  retrieve a UQLab model from the current session
%    myModel = UQ_GETMODEL returns the currently selected MODEL object
%    myModel from the UQLab session.
%
%    myModel = UQ_GETMODEL(MODELNAME) returns the MODEL object with the
%    specified name MODELNAME, if it exists in the UQLab session.
%    Otherwise, it returns an error.
%    
%    myModel = UQ_GETMODEL(N) returns the Nth MODEL object stored in the
%    UQLab session.
%
%    To print a list of the currently existing models, their corresponding
%    numbers and the currently selected one, use the <a href="matlab:help uq_listModels">uq_listModels</a> command.
%
%    See also: uq_createModel, uq_evalModel, uq_listModels, uq_selectModel,
%              uq_getInput, uq_getAnalysis 
%
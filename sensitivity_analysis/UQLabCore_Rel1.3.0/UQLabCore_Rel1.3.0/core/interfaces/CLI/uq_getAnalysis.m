% UQ_GETANALYSIS  retrieve a UQLab ANALYSIS object from the current session.
%    myAnalysis = UQ_GETANALYSIS returns the currently selected ANALYSIS object
%    myAnalysis from the UQLab session.
%
%    myAnalysis = UQ_GETANALYSIS(ANALYSISNAME) returns the ANALYSIS object with the
%    specified name ANALYSISNAME, if it exists in the UQLab session.
%    Otherwise, it returns an error.
%    
%    myAnalysis = UQ_GETANALYSIS(N) returns the Nth ANALYSIS object stored in the
%    UQLab session.
%
%    To print a list of the currently existing ANALYSIS objects, their 
%    numbers and the currently selected one, use the <a href="matlab:help uq_listAnalyses">uq_listAnalyses</a> command.
%
%    See also: uq_createAnalysis, uq_listAnalyses, uq_selectAnalysis, 
%              uq_getInput, uq_getModel
%
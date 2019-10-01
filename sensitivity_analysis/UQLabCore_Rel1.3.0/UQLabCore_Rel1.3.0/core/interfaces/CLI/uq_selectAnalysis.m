% UQ_SELECTANALYSIS   select an ANALYSIS object in the UQLab session.
%    UQ_SELECTANALYSIS interactively prompts the user to select one of the
%    available UQLab ANALYSIS objects stored in the current session. The
%    selected ANALYSIS is used by default by other UQLab commands, e.g.
%    <a href="matlab:help uq_getAnalysis">uq_getAnalysis</a>.
%
%    UQ_SELECTANALYSIS(ANALYSISNAME) selects the ANALYSIS object with property
%    'Name' equal to the specified ANALYSISNAME.
%
%    UQ_SELECTANALYSIS(N) selects the Nth created ANALYSIS.
%
%    myAnalysis = UQ_SELECTANALYSIS(...) also returns the selected ANALYSIS object
%    in the myAnalysis variable. 
%    
%    To print a list of the currently existing ANALYSIS objects, their 
%    numbers and the currently selected one, use the <a href="matlab:help uq_listAnalyses">uq_listAnalyses</a> command.
%
%    See also: uq_createAnalysis, uq_getAnalysis, uq_listAnalyses, 
%              uq_selectInput, uq_selectModel
%
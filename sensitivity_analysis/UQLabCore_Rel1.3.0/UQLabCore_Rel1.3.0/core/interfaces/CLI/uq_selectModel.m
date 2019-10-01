% UQ_SELECTMODEL   select a MODEL object in the UQLab session
%    UQ_SELECTMODEL interactively prompts the user to select one of the
%    available UQLab MODEL objects stored in the current session. The
%    selected MODEL is used by default by other UQLab commands, e.g.
%    <a href="matlab:help uq_evalModel">uq_evalModel</a> and <a href="matlab:help uq_getModel">uq_getModel</a>.
%
%    UQ_SELECTMODEL(MODELNAME) selects the MODEL object with property
%    'Name' equal to the specified MODELNAME.
%
%    UQ_SELECTMODEL(N) selects the Nth created MODEL.
%
%    myModel = UQ_SELECTMODEL(...) also returns the selected MODEL object
%    in the myModel variable. 
%    
%    To print a list of the currently existing MODEL objects, their 
%    numbers and the currently selected one, use the <a href="matlab:help uq_listModels">uq_listModels</a> command.
%
%    See also: uq_createModel, uq_evalModel, uq_getModel, uq_listModels,
%              uq_selectInput, uq_selectAnalysis 
%
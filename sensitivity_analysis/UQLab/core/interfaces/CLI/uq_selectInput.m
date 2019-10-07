% UQ_SELECTINPUT   select an INPUT object in the UQLab session.
%    UQ_SELECTINPUT   interactively prompts the user to select one of the
%    available UQLab INPUT objects stored in the current session. The
%    selected INPUT is used by default by other UQLab commands, e.g.
%    <a href="matlab:help uq_getSample">uq_getSample</a> and <a href="matlab:help uq_getInput">uq_getInput</a>.
%
%    UQ_SELECTINPUT(INPUTNAME) selects the INPUT object with property
%    'Name' equal to the specified INPUTNAME.
%
%    UQ_SELECTINPUT(N) selects the Nth created INPUT.
%
%    myInput = UQ_SELECTINPUT(...) also returns the selected INPUT object
%    in the myInput variable. 
%    
%    To print a list of the currently existing INPUT objects, their
%    numbers and the currently selected one, use the <a href="matlab:help uq_listInputs">uq_listInputs</a> command.
%
%    See also: uq_createInput, uq_getInput, uq_listInputs, uq_getSample,
%              uq_selectModel, uq_selectAnalysis 
%
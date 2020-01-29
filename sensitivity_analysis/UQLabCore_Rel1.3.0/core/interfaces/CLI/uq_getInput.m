% UQ_GETINPUT   retrieve a UQLab INPUT object from the current session.
%    myInput = UQ_GETINPUT returns the currently selected INPUT object
%    myInput from the UQLab session.
%
%    myInput = UQ_GETINPUT(INPUTNAME) returns the INPUT object with the
%    specified name INPUTNAME, if it exists in the UQLab session.
%    Otherwise, it returns an error.
%    
%    myInput = UQ_GETINPUT(N) returns the Nth INPUT object stored in the
%    UQLab session.
%
%    To print a list of the currently existing INPUT objects, their 
%    numbers and the currently selected one, use the <a href="matlab:help uq_listInputs">uq_listInputs</a> command.
%
%    See also: uq_createInput, uq_listInputs, uq_selectInput, uq_getSample,
%              uq_getModel, uq_getAnalysis 
%
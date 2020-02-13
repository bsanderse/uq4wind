% UQLAB   Initialize the UQLab uncertainty quantification software
%    UQLab is the uncertainty quantification software developed in the Matlab 
%    environment at ETH Zurich. It allows you to set up ingredients for a UQ 
%    analysis, namely to define a probabilistic INPUT model (input random 
%    variables), a computational MODEL, and to select the ANALYSIS you want 
%    to carry out, e.g. create a surrogate model, compute sensitivity indices, 
%    compute a probability of failure, etc.
%
%    Usage:
%    UQLAB - Initialize the UQLab framework and clear the current UQLab
%    session
%
%    UQLAB('SESSIONFILE.mat') - load the UQLab session file SESSIONFILE.mat
%    previously created with the uq_saveSession('SESSIONFILE.mat')command.
%    All the objects created can be accessed by using the <a href="matlab:help uq_getModel">uq_getModel</a>,
%    <a href="matlab:help uq_getInput">uq_getInput</a> and <a href="matlab:help uq_getAnslysis">uq_getAnalysis</a> commands.
%
%    See also: uq_saveSession, uq_createInput, uq_getInput, uq_createModel,
%              uq_getModel, uq_createAnalysis, uq_getAnalysis 
%
%    To access the list of available user manuals, please use the following 
%    command: <a href="matlab:uqlab -doc">uqlab -doc</a>
%
%    Copyright 2013-2017 ETH Zurich, all rights reserved
%

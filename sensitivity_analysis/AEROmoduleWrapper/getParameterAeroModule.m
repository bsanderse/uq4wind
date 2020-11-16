function [Parameters, UncertainInputs]  = getParameterAeroModule(turbineName)
% This routine is a wrapper that calls a Matlab function with the name of
% the turbine under consideration, as specified in the initialization file, 
% e.g. NM80_calibrate

% create a function handle of the turbine name,
% this function is supposed to exist in the case folder
% in this file all uncertain inputs should be defined
turbineData      = str2func(turbineName);
UncertainInputs  = turbineData();


% additionally, other (deterministic) parameters can be set, e.g. the
% turbine name
Parameters.turbineName = turbineName;

function Y = aero_module_calibration(X,P)
% Store the current working directory
workingDir = pwd;
%% Write random input data to aeromodule
writeAeroModuleInput(X,P)
%% Run the Aero module executable
chdir([pwd,'\AEROmodule\',P{29}]); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
chdir(workingDir)

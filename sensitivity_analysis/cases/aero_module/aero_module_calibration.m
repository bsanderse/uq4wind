function Y = aero_module_calibration(X,P)
% Store the current working directory
workingDir = pwd;
%% Write random input data to aeromodule
writeAeroModuleInput(X,P)
%% Run the Aero module executable
chdir([pwd,'\AEROmodule\',P{29}]); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
chdir(workingDir)
%% Read interpolated data
% This function writes the values in 'Y' handle for the QOI obtained
% from the Aero-Module run. 
data = read_interpolated_val(P{44}); 
if(strcmp(P{27},'Axial_Force'))
    Y = [data];
end
function Y = aero_module(X,P)
% Store the current working directory
workingDir = pwd; 

%% Write input data to aeromodule
writeAeroModuleInput(X,P)

%% Run the Aero module executable 
aeromodule_dir = fullfile(pwd,'AEROmodule',P{29},'current');
chdir(aeromodule_dir); % Go to the AERO module directory
disp('Running AeroModule...'); 
system('ECNAero.exe'); % Run the executable
chdir(workingDir)

%% Read output from the AERO module
% output directory of AeroModule
output_dir = fullfile(aeromodule_dir,'output');

% set-up handle to function which will return the output from the
% AeroModule at the appropriate location
outputFile = str2func(strcat(P{29},'_readoutput'));

% this calls for example NM80_calibrate_readoutput()
% or NM80_readoutput()
Y = outputFile(output_dir,P);

end
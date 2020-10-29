function Y = aero_module(X,P)
% Store the current working directory
workingDir = pwd; 

%% Write input data to aeromodule
writeAeroModuleInputReplacement(X,P)

%% Run the Aero module executable 
current_folder  = P.FixedParameters.current_folder;
chdir(current_folder); % Go to the AERO module directory
disp('Running AeroModule...'); 
system('ECNAero.exe'); % Run the executable
chdir(workingDir)

%% Read output from the AERO module
% output directory of AeroModule
output_dir = fullfile(current_folder,'output');

% set-up handle to function which will return the output from the
% AeroModule at the appropriate location
outputFile = str2func(strcat(P.FixedParameters.turbineName,'_readoutput'));

% this calls for example NM80_calibrate_readoutput()
% or NM80_readoutput()
Y = outputFile(output_dir,P);

end
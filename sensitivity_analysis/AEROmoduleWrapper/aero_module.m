function Y = aero_module(X,P)
% Store the current working directory
% workingDir = pwd; 

%% Write input data to aeromodule
writeAeroModuleInputReplacement(X,P)

%% Run the Aero module executable 
root_folder     = P.FixedParameters.root_folder;
current_folder  = P.FixedParameters.current_folder;
cur_dir         = fullfile(root_folder,current_folder);

chdir(fullfile(cur_dir)); % Go to the AERO module directory
disp('Running AeroModule...'); 
system('ECNAero.exe'); % Run the executable
chdir(root_folder)

%% Read output from the AERO module
% output directory of AeroModule
output_dir = fullfile(cur_dir,'output');

% set-up handle to function which will return the output from the
% AeroModule at the appropriate location
outputFile = str2func(strcat(P.FixedParameters.turbineName,'_readoutput'));

% this calls for example NM80_calibrate_readoutput()
% or NM80_readoutput()
Y = outputFile(output_dir,P);

end
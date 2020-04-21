function Y = aero_module_axial(X,P)
% Store the current working directory
workingDir = pwd; 
%% Write random input data to aeromodule
writeAeroModuleInput(X,P)
%% Run the Aero module executable                            
chdir([pwd,'\AEROmodule\',P{29}]); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
chdir(workingDir)
%% Read output from the AERO module
output_axial = readtable(pwd,'\AEROmodule\',P{29},'\output\output.dat'); % location of the AeroPower.dat output file
output = importfile(output_axial, 2);
% return output 
if(strcmp(P{27},'rad'))
    Y = output.r_i;
elseif(strcmp(P{27},'axial'))
    Y = output.aero_data;  
end
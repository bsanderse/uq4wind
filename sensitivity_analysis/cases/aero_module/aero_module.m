function Y = aero_module(X,P)
% Store the current working directory
workingDir = pwd; 
%% Write random input data to aeromodule
writeAeroModuleInput(X,P)
%% Run the Aero module executable                            
chdir([pwd,'\AEROmodule\']); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
chdir(workingDir)
%% Read output from the AERO module
filename = [pwd,'\AEROmodule\','output\AeroPower.dat']; % location of the AeroPower.dat output file
[Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});
% return output 
if(strcmp(P{27},'Power'))
    Y = mean(PowerWatt);
elseif(strcmp(P{27},'Axial_Force'))
    Y = mean(Axial_ForceN);  
end
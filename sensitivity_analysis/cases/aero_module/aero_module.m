function Y = aero_module(X,P)
% Store the current working directory
workingDir = pwd; 

%% Write input data to aeromodule
writeAeroModuleInput(X,P)

%% Run the Aero module executable                            
chdir([pwd,'\AEROmodule\',P{29}]); % Go to the AERO module directory
disp('Running AeroModule...');
test = system('ECNAero.exe'); % Run the executable
chdir(workingDir)

%% Read output from the AERO module
% output directory of AeroModule
output_dir = strcat(pwd,'\AEROmodule\',P{29},'\output\');

% return output depending on the QoI we are interested in
switch P{27}
    case 'Power'
        filename = strcat(output_dir,'AeroPower.dat'); % location of the AeroPower.dat output file
        [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});
        Y = mean(PowerWatt);
    case 'Axial_Force'
        filename = strcat(output_dir,'AeroPower.dat'); % location of the AeroPower.dat output file
        [Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});
        Y = mean(Axial_ForceN); 
    case 'Sectional_normal_force'   
        Y = read_interpolated_val(strcat(output_dir,'B1n_BEM.txt')); 
    otherwise
        error(strcat('QoI type unknown; check the turbine file ',P{29}));
end
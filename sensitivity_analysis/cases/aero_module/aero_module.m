function Y = aero_module(X,P)

% Store the current working directory
workingDir = pwd; 

% Location where the AERO modules 
folder ='C:\Users\pkumar\Dropbox\WindTrue\ECNAero2CWI\';

% Get the random samples of twist and chord
twist = computeTwist(1, randVec, pc, 0); % computeTwist routine uses the specifications of NM80 turbine by default
chord = computeChord(1, randVec, pc, 0); % computeChord routine uses the specifications of NM80 turbine by default

filename = [folder,'output\AeroPower.dat']; % location of the AeroPower.dat output file
% Loop through the samples

% Write the random data into the input file of AERO module
writeAeroModuleInput(P{1},P{2},P{10}, P{3}, chord, ...
                     P{5}, twist, P{7}, P{8}, P{9},... 
                     P{11}, P{12}, P{13}, P{14}, P{17}, ...
                     P{15}, P{19}, P{16}, P{18}, P{20},...
                     P{21}, P{22}, folder)
                           
chdir(folder); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
% read output from the AERO module
[Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});


chdir(workingDir)

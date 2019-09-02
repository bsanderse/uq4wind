function Y = aero_module(X,P)
% Store the current working directory
workingDir = pwd; 

%% Get one random samples for twist and chord
randVecTwist = zeros(1,P{26}); 
randVecTwist(P{30}) = X(1:length(P{30}));
randVecChord = zeros(1,P{27});
randVecChord(P{31}) = X(length(P{30})+1:end-1);
twist = computeTwist(1, randVecTwist, P{28}, 0); % computeTwist routine uses the specifications of NM80 turbine by default
chord = computeChord(1, randVecChord, P{29}, 0); % computeChord routine uses the specifications of NM80 turbine by default

%% Get one random samples for YAWANGLE
YAWANGLE = P{20} + P{32}.*X(end);
%% Write the random data into the input file of AERO module
writeAeroModuleInput(P{1},  P{2},  P{10}, P{3},  chord, ...
                     P{5},  twist, P{7},  P{8},  P{9}, ... 
                     P{11}, P{12}, P{13}, P{14}, P{17}, ...
                     P{15}, P{19}, P{16}, P{18}, YAWANGLE, ...
                     P{21}, P{22}, P{25})
%% Run the Aero module executable                            
chdir(P{25}); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
%% Read output from the AERO module
filename = [P{25},'output\AeroPower.dat']; % location of the AeroPower.dat output file
[Times,Azimuthdeg,PowerWatt,Axial_ForceN] = AeroPower(filename, P{23}, P{24});
chdir(workingDir)
% return output 
Y = mean(PowerWatt);
% Y = Axial_ForceN;  
function Y = aero_module(X,P)
% Store the current working directory
workingDir = pwd; 

% Location where the AERO modules 
folder ='C:\Users\pkumar\Dropbox\WindTrue\ECNAero2CWI\';
% Get the random samples of twist and chord

twist = computeTwist(samples, 1); % computeTwist routine uses the specifications of NM80 turbine by default
chord = computeChord(samples, 1); % computeChord routine uses the specifications of NM80 turbine by default

filename = [folder,'output\AeroPower.dat']; % location of the AeroPower.dat output file

% Loop through the samples
for i = 1:samples
% Write the random data into the input file of AERO module
writeAeroModuleInput(AEROMODEL,TURBINETYPE,vectorLength,...
                              zB, chord(i,:), t_by_c, twist(i,:), C14, xB, yB,... 
                              BLADELENGTH, BLADEROOT, HUBHEIGHT, TILTANGLE, RPM, ...
                              PITCHANGLE, TIMESTEP, XNAC2HUB, TEND, YAWANGLE,...
                              NROFBEMELEMENTS, ZNAC2HUB, folder)
                           
chdir(folder); % Go to the AERO module directory
system('ECNAero.exe') % Run the executable
% read output from the AERO module
[Times,Azimuthdeg,PowerWatt(i,:),Axial_ForceN(i,:)] = AeroPower(filename, startRow, endRow);
end

% Also compute the baseline case
writeAeroModuleInput(AEROMODEL,TURBINETYPE,vectorLength,...
                              zB, ref_chord, t_by_c, ref_twist, C14, xB, yB,... 
                              BLADELENGTH, BLADEROOT, HUBHEIGHT, TILTANGLE, RPM, ...
                              PITCHANGLE, TIMESTEP, XNAC2HUB, TEND, YAWANGLE,...
                              NROFBEMELEMENTS, ZNAC2HUB, folder)
                           
system('ECNAero.exe') % Run the executable
[Times,Azimuthdeg,ref_PowerWatt,ref_Axial_ForceN] = AeroPower(filename, startRow, endRow);

figure % Plot the samples of Power
plot(Times,PowerWatt','linewidth',1,'color','k','HandleVisibility','off')
hold on
plot(Times,ref_PowerWatt','linewidth',2)
title('Samples of power (Watt)')
legend('Baseline power')
xlabel('Times (sec)')
ylabel('Power (Watt)')

figure % Plot the samples of Axial force
plot(Times,Axial_ForceN','linewidth',1,'color','k','HandleVisibility','off')
hold on
plot(Times,ref_Axial_ForceN','linewidth',2)
title('Samples of axial force (N)')
xlabel('Times (sec)')
ylabel('Axial force (N)')
legend('Baseline axial force')
chdir(workingDir)

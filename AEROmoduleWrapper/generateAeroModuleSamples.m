% This routine compute samples of power output and axial force with uncertain chord and twist inputs
% using the ECN AERO module software package 

samples = 10; % Number of samples 
folder ='C:\Users\pkumar\Dropbox\WindTrue\ECNAero2CWI\'; % Location where the AERO modules  

% Variables of input file extracted from reference test case from DANAERO
% turbuine NM80
AEROMODEL = 1;
TURBINETYPE = 1;

zB = [0 2 4 6 8 10 12 14 16 18 20 22 ...
      24 26 28 30 32 34 36 37 38 38.4 38.8];

chord = [2.42 2.48 2.65 2.81 2.98 3.14 3.17 2.99 2.79 2.58 2.38 ...
         2.21 2.06 1.92 1.8 1.68 1.55 1.41 1.18 0.98 0.62 0.48 0.07]; % Replace with random chord vector using NURBS

t_by_c = [0.9999 0.9641 0.8053 0.6508 0.5167 0.403 0.3253 0.284 0.2562 0.2377 0.2225 ...
          0.2099 0.2003 0.194 0.1903 0.1879 0.186 0.1839 0.1795 0.1739 0.1633 0.157 0.1484]; % Replace with ?? 

twist = [0 5.37 6.69 7.9 9.11 10.19 9.39 7.16 5.45 4.34 3.5 2.86 ... 
         2.31 1.77 1.28 0.9 0.55 0.23 0.03 0.02 0.93 2.32 6.13]; % Replace with random twist vector using NURBS

C14 = [-25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 ...
       -25 -25 -25 -25 -25 -25 -25 -25 -25];

xB = [0 0.00145 0.00614 0.01254 0.02073 0.0169 -0.0023 -0.0376614 -0.07744 -0.12355 -0.17826 ...
     -0.24389 -0.32359 -0.41925 -0.53316 -0.66737 -0.82386 -1.00305 -1.20272 -1.3091 -1.41608 -1.45947 -1.5023];

yB = [0 0.01545 0.05235 0.0896 0.12624 0.16279 0.19638 0.20193 0.2072 0.21224 0.21281 0.20968 0.20203 ...
      0.19167 0.17959 0.16755 0.15536 0.14092 0.11843 0.09797 0.06222 0.04084 0.00036];

vectorLength = length(zB);
BLADELENGTH = 38.8;
BLADEROOT = 1.24;
HUBHEIGHT = 57.19;
TILTANGLE = 0.0;
PITCHANGLE = 0.15;
XNAC2HUB = -4.03;
RPM = 12.3;
TEND = 60;
TIMESTEP = 0.135501355;

%[TIMESTEP,TEND] = wakepoints(RPM); 
YAWANGLE = 0.0;
NROFBEMELEMENTS = 26;
ZNAC2HUB = 1.6;

% Write the random data into the input file of AERO module
writeAeroModuleInput(AEROMODEL,TURBINETYPE,vectorLength,...
                              zB, chord, t_by_c, twist, C14, xB, yB,... 
                              BLADELENGTH, BLADEROOT, HUBHEIGHT, TILTANGLE, RPM, ...
                              PITCHANGLE, TIMESTEP, XNAC2HUB, TEND, YAWANGLE,...
                              NROFBEMELEMENTS, ZNAC2HUB, folder)



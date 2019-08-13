function P = getParameterAeroModule()
% This routine returns the deterministic parameters set 'P' for the input
% file of AERO module software.

% Define parameter set as empty cell for the uq_lab
P ={}; 

%% Variables of input file extracted from reference test case from DANAERO
% turbuine NM80
AEROMODEL = 1;
TURBINETYPE = 1;

zB = [0 2 4 6 8 10 12 14 16 18 20 22 ...
      24 26 28 30 32 34 36 37 38 38.4 38.8];
  
ref_chord = [2.42 2.48 2.65 2.81 2.98 3.14 3.17 2.99 2.79 2.58 2.38 ...
         2.21 2.06 1.92 1.8 1.68 1.55 1.41 1.18 0.98 0.62 0.48 0.07]; % Baseline chord, not needed in the simulation...
     
t_by_c = [0.9999 0.9641 0.8053 0.6508 0.5167 0.403 0.3253 0.284 0.2562 0.2377 0.2225 ...
          0.2099 0.2003 0.194 0.1903 0.1879 0.186 0.1839 0.1795 0.1739 0.1633 0.157 0.1484]; % Thickness by chord ratio should also be uncertain, currently considered deterministic
      
ref_twist = [0 5.37 6.69 7.9 9.11 10.19 9.39 7.16 5.45 4.34 3.5 2.86 ... 
         2.31 1.77 1.28 0.9 0.55 0.23 0.03 0.02 0.93 2.32 6.13]; % Baseline twist, not needed in the simulation...
     
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
TEND = 30;
TIMESTEP = 0.135501355;
%[TIMESTEP,TEND] = wakepoints(RPM); % Routine to compute the TIMESTEP and  TEND using RPM
YAWANGLE = 0.0;
NROFBEMELEMENTS = 26;
ZNAC2HUB = 1.6;

%% Other parameters
% Top and the end of AeroPower.dat file
startRow = 2;
endRow = round(TEND/TIMESTEP)+ startRow; % depends on the TEND and TIMESTEP
% Location where the AERO module is stored
folder ='C:\Users\pkumar\Dropbox\WindTrue\ECNAero2CWI\';

%% Define properties of uncerain input
% Number of control points twist and chord. Set heuristically.
NRCP_TWIST = 10;
NRCP_CHORD = 9; 

% fraction of perturbation for each control points, 0.1 corresponds to plus
% minus 5% perturbation on the baseline values
PERTURBATION_TWIST = 0.1*ones(1,NRCP_TWIST);  
PERTURBATION_CHORD = 0.1*ones(1,NRCP_CHORD);  

% Index at which we want to introduce the uncertainty. This is to control
% the number of uncertain paramters. We only introduce uncertainties in the
% "important" control points.
INDEX_TWIST = [2 3 4]; % INDEX_TWIST = 1:NRCP_TWIST
INDEX_CHORD = [2 3 4]; % INDEX_CHORD = 1:NRCP_CHORD

%% Populate the elements of P. The order of elements in P should not be changed
P{1} = AEROMODEL; 
P{2} = TURBINETYPE;
P{3} = zB;
P{4} = ref_chord;
P{5} = t_by_c;
P{6} = ref_twist;
P{7} = C14;
P{8} = xB;
P{9} = yB;
P{10} = vectorLength;
P{11} = BLADELENGTH;
P{12} = BLADEROOT;
P{13} = HUBHEIGHT;
P{14} = TILTANGLE;
P{15} = PITCHANGLE;
P{16} = XNAC2HUB;
P{17} = RPM;
P{18} = TEND;
P{19} = TIMESTEP;
P{20} = YAWANGLE;
P{21} = NROFBEMELEMENTS;
P{22} = ZNAC2HUB;
P{23} = startRow;
P{24} = endRow;
P{25} = folder;
P{26} = NRCP_TWIST;
P{27} = NRCP_CHORD;
P{28} = PERTURBATION_TWIST;
P{29} = PERTURBATION_CHORD;
P{30} = INDEX_TWIST;
P{31} = INDEX_CHORD;


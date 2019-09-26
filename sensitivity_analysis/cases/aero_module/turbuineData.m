function [AEROMODEL,TURBINETYPE,zB, ref_chord, t_by_c,ref_twist, C14, xB, yB, vectorLength, ...
          BLADELENGTH, BLADEROOT, HUBHEIGHT, TILTANGLE, PITCHANGLE, XNAC2HUB, ...
          RPM, TEND, TIMESTEP, YAWANGLE, NROFBEMELEMENTS, ZNAC2HUB, Input, uncertain_params, QoI]  = turbuineData()
%% Variables of input file extracted from reference test case from DANAERO turbuine NM80
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

%% Define properties of uncertain input in the UQLab format. 
% We define this for all possible uncertain inputs and finally in the 
% variable "uncertain_params" we specify which variables to be considered                                    
% for the sensitivity analysis. NOTE that the ordering of variables is not important.

counter = 0; % To keep track of the number of uncertain variables. The value should be increase by one after adding every new random parameter. 

%% ===============TWIST====================
% Using 12 control points for twist chosen heuristically. Variable "Input.Marginals.Index" stores the index of
% control points. Here we assume same marginal for all control points.
NRCP_TWIST = 12;
for i=1:NRCP_TWIST
    counter = counter+1;
    Input.Marginals(counter).Name = 'Twist';
    Input.Marginals(counter).Index = i;
    Input.Marginals(counter).Type = 'Uniform'; 
    Input.Marginals(counter).Parameters = [-0.5 0.5];
    Input.Marginals(counter).Bounds = [-0.5 0.5];
end

%% ===============CHORD====================
% Using 9 control points for chord chosen heuristically. Variable "Input.Marginals.Index" stores the index of
% control points.
NRCP_CHORD = 9;
for i=1:NRCP_CHORD
    counter = counter+1;
    Input.Marginals(counter).Name = 'Chord';
    Input.Marginals(counter).Index = i;
    Input.Marginals(counter).Type = 'Uniform'; 
    Input.Marginals(counter).Parameters = [-0.5 0.5];
    Input.Marginals(counter).Bounds = [-0.5 0.5];
end

%% ===============THICKNESS====================
% Using 9 control points for thickness chosen heuristically. Variable "Input.Marginals.Index" stores the index of
% control points.
NRCP_THICKNESS = 9;
for i=1:NRCP_THICKNESS
    counter = counter+1;
    Input.Marginals(counter).Name = 'Thickness';
    Input.Marginals(counter).Index = i;
    Input.Marginals(counter).Type = 'Uniform'; 
    Input.Marginals(counter).Parameters = [-0.5 0.5];
    Input.Marginals(counter).Bounds = [-0.5 0.5];
end

%%=======================YAW====================
% Truncated Gaussian
YAW_Std = 2;  % Standard deviation
YAW_LB = -10; % Lower bound of trucated Gaussian distribution
YAW_UB = 10;  % Upper bound of trucated Gaussian distribution
counter = counter+1;
Input.Marginals(counter).Name = 'YAW';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [YAWANGLE, YAW_Std];
Input.Marginals(counter).Bounds = [YAW_LB YAW_UB]; 

%%=======================PITCH==================
%%=======================RPM====================
%%=======================WINDSPEED==============

% Specify uncertain parameters {name,index}to consider in the sensitivity analysis
uncertain_params = {{'Twist',2},{'Twist',3},{'Twist',4},{'Twist',5},{'Twist',6},{'Twist',7},{'Chord',2},{'Chord',4},{'Chord',6},{'Chord',8}, ...
                    {'Thickness',2},{'Thickness',3},{'Thickness',4},{'Thickness',5},{'YAW',0}};

% Specify quantity of interest
QoI = 'Power'; % 'Axial_Force'




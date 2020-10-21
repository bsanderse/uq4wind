function [ Input, uncertain_params, QoI, Inp]  = NM80()

%% Specify uncertain parameters to be considered in the sensitivity analysis
% The parameter should be defined in the following format {name,index,rel_perturbation} where
% rel_pertubation defines the amount of relative perturbation for B-spline
% curves. For example, {'Twist',2,0.2} defines the uncertainty in the
% second control point of Twist curve and 0.2 means a relative uncertainty of plus minus 10% for this  control point. 
% This parameter may not be required for other scalar random
% variables. 

% NOTE: the index associated to Twist, Thickness and Chord is the
% associated to the zB vector
% the index associated to Cl, Cd, and Cm is associated to the available airfoil
% polars, which is generally at a different section

% uncertain_params = {{'Twist',2,0.2},{'Twist',3,0.2},{'Twist',4,0.2},{'Twist',5,0.2},{'Twist',6,0.2},{'Twist',7,0.2},...
%                    {'Chord',2,0.2},{'Chord',4,0.2},{'Chord',6,0.2},{'Chord',8,0.2}, ...
%                    {'Thickness',2,0.2},{'Thickness',3,0.2},{'Thickness',4,0.2},{'Thickness',5,0.2},...
%                    {'YAW','',''},{'WINDSPEED','',''},{'RPM','',''},{'PITCHANGLE','',''}, ...
%                    {'CL',1, 0.2}, {'CL',2, 0.2}, {'CL',3, 0.2},{'CL',4,0.2}, ...
%                    {'CD',1, 0.2}, {'CD',2, 0.2}, {'CD',3, 0.2},{'CD',4,0.2}, ...
%                    {'CM',1, 0.2}, {'CM',2, 0.2}, {'CM',3, 0.2},{'CM',4,0.2}
%                    {'DYNSTALLTYPE','',''}, {'CORR3DTYPE','',''},...
%                    {'BL_A1','',''},{'BL_A2','',''},{'BL_b1','',''},{'BL_b2','',''},...
%                    {'BL_Ka','',''},{'BL_Tp','',''},{'BL_Tf','',''},{'BL_Tv','',''},{'BL_Tvl','',''},{'BL_Acd','',''}};

 uncertain_params = {{'CL',1, 0.2}, {'CL',2, 0.2}, {'CL',3, 0.2}, {'CL',4,0.2}};

%uncertain_params = {{'BL_A1','',''},{'BL_A2','',''},{'BL_b1','',''},{'BL_b2','',''},...
%    {'BL_Ka','',''},{'BL_Tp','',''},{'BL_Tf','',''},{'BL_Tv','',''}};
                
% Specify quantity of interest
QoI = 'Power'; % 'Axial_Force' or 'Power'


%% Variables of input file extracted from reference test case from DANAERO turbine NM80
Inp.AEROMODEL = 1;
Inp.TURBINETYPE = 1;

Inp.zB = [0 2 4 6 8 10 12 14 16 18 20 22 ...
      24 26 28 30 32 34 36 37 38 38.4 38.8];
  
Inp.ref_chord = [2.42 2.48 2.65 2.81 2.98 3.14 3.17 2.99 2.79 2.58 2.38 ...
         2.21 2.06 1.92 1.8 1.68 1.55 1.41 1.18 0.98 0.62 0.48 0.07]; % Baseline chord, not needed in the simulation...
     
Inp.t_by_c = [0.9999 0.9641 0.8053 0.6508 0.5167 0.403 0.3253 0.284 0.2562 0.2377 0.2225 ...
          0.2099 0.2003 0.194 0.1903 0.1879 0.186 0.1839 0.1795 0.1739 0.1633 0.157 0.1484]; % Thickness by chord ratio should also be uncertain, currently considered deterministic
      
Inp.ref_twist = [0 5.37 6.69 7.9 9.11 10.19 9.39 7.16 5.45 4.34 3.5 2.86 ... 
         2.31 1.77 1.28 0.9 0.55 0.23 0.03 0.02 0.93 2.32 6.13]; % Baseline twist, not needed in the simulation...
     
Inp.C14 = [-25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 -25 ...
       -25 -25 -25 -25 -25 -25 -25 -25 -25];
   
Inp.xB = [0 0.00145 0.00614 0.01254 0.02073 0.0169 -0.0023 -0.0376614 -0.07744 -0.12355 -0.17826 ...
     -0.24389 -0.32359 -0.41925 -0.53316 -0.66737 -0.82386 -1.00305 -1.20272 -1.3091 -1.41608 -1.45947 -1.5023];

Inp.yB = [0 0.01545 0.05235 0.0896 0.12624 0.16279 0.19638 0.20193 0.2072 0.21224 0.21281 0.20968 0.20203 ...
      0.19167 0.17959 0.16755 0.15536 0.14092 0.11843 0.09797 0.06222 0.04084 0.00036];

Inp.vectorLength = length(Inp.zB);
Inp.BLADELENGTH = 38.8;
Inp.BLADEROOT = 1.24;
Inp.HUBHEIGHT = 57.19;
Inp.TILTANGLE = 0.0;
Inp.PITCHANGLE = 0.15;
Inp.XNAC2HUB = -4.03;
Inp.RPM = 12.3;
Inp.TBEGIN = 0.0;
Inp.TEND = 30;
Inp.TIMESTEP = 0.135501355;
%[TIMESTEP,TEND] = wakepoints(RPM); % Routine to compute the TIMESTEP and  TEND using RPM
Inp.YAWANGLE = 0.0;
Inp.NROFBEMELEMENTS = 26;
Inp.ZNAC2HUB = 1.6;
Inp.WINDSPEED = 6.1;
Inp.DYNSTALLTYPE = 3; % 0: No DS 1:Snel1 2: Snel2 3:B-Leishmann 4: Onera
Inp.CORR3DTYPE = 0;
Inp.AIRDENSITY = 1.231;
Inp.TOWERBASERADIUS = 0.0001;
Inp.TOWERTOPRADIUS = 0.0001;
% BL model parameters
Inp.BL_A1 = 0.30;
Inp.BL_A2 = 0.70;
Inp.BL_b1 = 0.14;
Inp.BL_b2 = 0.53;
Inp.BL_Ka = 0.75;
Inp.BL_Tp = 1.50;
Inp.BL_Tf = 5.00;
Inp.BL_Tv = 6.00;
Inp.BL_Tvl = 5.00;
Inp.BL_Acd = 0.13;

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

%% =======================YAW====================
% Truncated Gaussian
YAWANGLE   = Inp.YAWANGLE;
YAW_Std = 2;  % Standard deviation
YAW_LB = -10; % Lower bound of trucated Gaussian distribution
YAW_UB = 10;  % Upper bound of trucated Gaussian distribution
counter = counter+1;
Input.Marginals(counter).Name = 'YAW';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [YAWANGLE, YAW_Std];
Input.Marginals(counter).Bounds = [YAW_LB YAW_UB]; 

%% =======================PITCH==================
% Truncated Gaussian
PITCHANGLE   = Inp.PITCHANGLE;
PITCHANGLE_Std = 1;  % Standard deviation
PITCHANGLE_LB = -2; % Lower bound of trucated Gaussian distribution
PITCHANGLE_UB = 2;  % Upper bound of trucated Gaussian distribution
counter = counter+1;
Input.Marginals(counter).Name = 'PITCHANGLE'; 
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [PITCHANGLE, PITCHANGLE_Std];
Input.Marginals(counter).Bounds = [PITCHANGLE_LB PITCHANGLE_UB];

%% =======================RPM====================
% Truncated Gaussian
RPM = Inp.RPM;
RPM_Std = 1;  % Standard deviation
RPM_LB = 10; % Lower bound of trucated Gaussian distribution
RPM_UB = 14;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'RPM';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [RPM, RPM_Std];
Input.Marginals(counter).Bounds = [RPM_LB RPM_UB]; 

%% =======================WINDSPEED==============
WindSpeed_scale = 6.1;
WindSpeed_shape = 50;
counter = counter + 1;
Input.Marginals(counter).Name = 'WINDSPEED';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Weibull'; 
Input.Marginals(counter).Parameters = [WindSpeed_scale WindSpeed_shape]; % scale and shape parameter
Input.Marginals(counter).Bounds = ''; % No bound needed for Weibull 

%% =======================DYNSTALLTYPE==============
% Discrete variable with values 0,1,2,3,4 with 0: No DS 1:Snel1 2: Snel2
% 3:B-Leishmann 4: Onera. We use uniform random variable from [0 5] sample models with equal probability

counter = counter+1;
Input.Marginals(counter).Name = 'DYNSTALLTYPE';
Input.Marginals(counter).Index = '';
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [0 20-10^-20]; % subtract a small bias as floor operator should result in 0,1,2,3,4 
Input.Marginals(counter).Bounds = [0 20-10^-20];


%% =======================CORR3DTYPE==============
% Discrete variable with values 0,1  with 0: No correction 1: Snel. 
% We use uniform random variable from [0 2] sample models with equal probability
counter = counter+1;
Input.Marginals(counter).Name = 'CORR3DTYPE';
Input.Marginals(counter).Index = '';
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [0 8-10^-20]; % subtract a small bias as floor operator should result in 0,1
Input.Marginals(counter).Bounds = [0 8-10^-20];


%% ====================Polars=====================
% Import the reference polar curves
% using fullfile works both on Windows and Linux and Mac
section03_file = fullfile('..','..','AEROmodule','NM80_calibrate','reference','section03_ref.dat');
section05_file = fullfile('..','..','AEROmodule','NM80_calibrate','reference','section05_ref.dat');
section08_file = fullfile('..','..','AEROmodule','NM80_calibrate','reference','section08_ref.dat');
section10_file = fullfile('..','..','AEROmodule','NM80_calibrate','reference','section10_ref.dat');
[alpha03, CL_section03, CD_section03, CM_section03] = importPolars(section03_file);
[alpha05, CL_section05, CD_section05, CM_section05] = importPolars(section05_file);
[alpha08, CL_section08, CD_section08, CM_section08] = importPolars(section08_file);
[alpha10, CL_section10, CD_section10, CM_section10] = importPolars(section10_file);
alpha03_pert = find(alpha03<=50 & alpha03>=-10);
alpha05_pert = find(alpha05<=50 & alpha05>=-10);
alpha08_pert = find(alpha08<=50 & alpha08>=-10);
alpha10_pert = find(alpha10<=50 & alpha10>=-10);
% Important: the format of POLARS cell defined below should not be changed 
Inp.POLARS = {4,... % Number of polar files
    {'section03.dat', 'section05.dat','section08.dat', 'section10.dat'}, ... % Name of polar files
    {'Section03', 'Section05', 'Section08', 'Section10'}, ... % Airfoil_Name
    {0.333, 0.243, 0.197, 0.187}, ... % thickness by chord ratio
    {1.0E+07, 1.0E+07, 1.0E+07, 1.0E+07,}, ... % Reynolds numbers
    {CL_section03, CD_section03, CM_section03, alpha03, alpha03_pert}, ... % Cl, Cd, Cm data for section03
    {CL_section05, CD_section05, CM_section05, alpha05, alpha05_pert}, ...
    {CL_section08, CD_section08, CM_section08, alpha08, alpha08_pert}, ...
    {CL_section10, CD_section10, CM_section10, alpha10, alpha10_pert}, ...
    };

% Define PDF for CL
counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Index = 1; % Corresponds to section 3
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Index = 2; % Corresponds to section 5
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Index = 3; % Corresponds to section 8
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Index = 4; % Corresponds to section 10
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

% Define PDF for CD
counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Index = 1; % Corresponds to section 3
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Index = 2; % Corresponds to section 5
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Index = 3; % Corresponds to section 8
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Index = 4; % Corresponds to section 10
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

% Define PDF for CM
counter = counter+1;
Input.Marginals(counter).Name = 'CM';
Input.Marginals(counter).Index = 1; % Corresponds to section 3
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CM';
Input.Marginals(counter).Index = 2; % Corresponds to section 5
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CM';
Input.Marginals(counter).Index = 3; % Corresponds to section 8
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CM';
Input.Marginals(counter).Index = 4; % Corresponds to section 10
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-0.5 0.5];
Input.Marginals(counter).Bounds = [-0.5 0.5];


%% Inputs for BL model parameters
% =======================BL_A1====================
% Truncated Gaussian
BL_A1   = Inp.BL_A1;
BL_A1_Std = 1;  % Standard deviation
BL_A1_LB = BL_A1 - BL_A1*0.1; % Lower bound of trucated Gaussian distribution
BL_A1_UB = BL_A1 + BL_A1*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_A1';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_A1, BL_A1_Std];
Input.Marginals(counter).Bounds = [BL_A1_LB BL_A1_UB]; 

% =======================BL_A2====================
% Truncated Gaussian
BL_A2   = Inp.BL_A2;
BL_A2_Std = 1;  % Standard deviation
BL_A2_LB = BL_A2 - BL_A2*0.1; % Lower bound of trucated Gaussian distribution
BL_A2_UB = BL_A2 + BL_A2*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_A2';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_A2, BL_A2_Std];
Input.Marginals(counter).Bounds = [BL_A2_LB BL_A2_UB]; 

% =======================BL_b1====================
% Truncated Gaussian
BL_b1   = Inp.BL_b1;
BL_b1_Std = 1;  % Standard deviation
BL_b1_LB = BL_b1 - BL_b1*0.1; % Lower bound of trucated Gaussian distribution
BL_b1_UB = BL_b1 + BL_b1*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_b1';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_b1, BL_b1_Std];
Input.Marginals(counter).Bounds = [BL_b1_LB BL_b1_UB];

% =======================BL_b2====================
% Truncated Gaussian
BL_b2   = Inp.BL_b2;
BL_b2_Std = 1;  % Standard deviation
BL_b2_LB = BL_b2 - BL_b2*0.1; % Lower bound of trucated Gaussian distribution
BL_b2_UB = BL_b2 + BL_b2*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_b2';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_b2, BL_b2_Std];
Input.Marginals(counter).Bounds = [BL_b2_LB BL_b2_UB];

% =======================BL_Ka====================
% Truncated Gaussian
BL_Ka   = Inp.BL_Ka;
BL_Ka_Std = 1;  % Standard deviation
BL_Ka_LB = BL_Ka - BL_Ka*0.1; % Lower bound of trucated Gaussian distribution
BL_Ka_UB = BL_Ka + BL_Ka*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Ka';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Ka, BL_Ka_Std];
Input.Marginals(counter).Bounds = [BL_Ka_LB BL_Ka_UB];

% =======================BL_Tp====================
% Truncated Gaussian
BL_Tp   = Inp.BL_Tp;
BL_Tp_Std = 1;  % Standard deviation
BL_Tp_LB = BL_Tp - BL_Tp*0.1; % Lower bound of trucated Gaussian distribution
BL_Tp_UB = BL_Tp + BL_Tp*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Tp';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Tp, BL_Tp_Std];
Input.Marginals(counter).Bounds = [BL_Tp_LB BL_Tp_UB];

 
% =======================BL_Tf====================
% Truncated Gaussian
BL_Tf   = Inp.BL_Tf;
BL_Tf_Std = 1;  % Standard deviation
BL_Tf_LB = BL_Tf - BL_Tf*0.1; % Lower bound of trucated Gaussian distribution
BL_Tf_UB = BL_Tf + BL_Tf*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Tf';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Tf, BL_Tf_Std];
Input.Marginals(counter).Bounds = [BL_Tf_LB BL_Tf_UB];

% =======================BL_Tv====================
% Truncated Gaussian
BL_Tv   = Inp.BL_Tv;
BL_Tv_Std = 1;  % Standard deviation
BL_Tv_LB = BL_Tv - BL_Tv*0.1; % Lower bound of trucated Gaussian distribution
BL_Tv_UB = BL_Tv + BL_Tv*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Tv';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Tv, BL_Tv_Std];
Input.Marginals(counter).Bounds = [BL_Tv_LB BL_Tv_UB];

% =======================BL_Tvl====================
% Truncated Gaussian
BL_Tvl   = Inp.BL_Tvl;
BL_Tvl_Std = 1;  % Standard deviation
BL_Tvl_LB = BL_Tvl - BL_Tvl*0.1; % Lower bound of trucated Gaussian distribution
BL_Tvl_UB = BL_Tvl + BL_Tvl*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Tv1';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Tvl, BL_Tvl_Std];
Input.Marginals(counter).Bounds = [BL_Tvl_LB BL_Tvl_UB];

% =======================BL_Acd====================
% Truncated Gaussian
BL_Acd   = Inp.BL_Acd;
BL_Acd_Std = 1;  % Standard deviation
BL_Acd_LB = BL_Acd - BL_Acd*0.1; % Lower bound of trucated Gaussian distribution
BL_Acd_UB = BL_Acd + BL_Acd*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Acd';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Acd, BL_Acd_Std];
Input.Marginals(counter).Bounds = [BL_Acd_LB BL_Acd_UB];




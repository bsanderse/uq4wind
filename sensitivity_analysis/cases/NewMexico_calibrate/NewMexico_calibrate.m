function [Input, QoI]  = NewMexico_calibrate()

QoI = 'Sectional_normal_force';

counter = 1;

RPM     = 425.1;
RPM_Std = 1;  % Standard deviation
RPM_LB = 10; % Lower bound of trucated Gaussian distribution
RPM_UB = 14;  % Upper bound of trucated Gaussian distribution 
Input.Marginals(counter).Name = 'RPM';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [RPM, RPM_Std];
Input.Marginals(counter).Bounds = [RPM_LB RPM_UB]; 

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'risoA121.dat';
Input.Marginals(counter).AlphaPert = [-10 30];
Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'du91w250.dat';
Input.Marginals(counter).AlphaPert = [-10 30];
Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
Input.Marginals(counter).AlphaPert = [-10 30];
Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CD';
% Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
% Input.Marginals(counter).AlphaPert = [0 20];
% Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];


end


%% ===============Variables for user selection====================

% uncertain_params = {{'CL',1, 0.2}, {'CL',3,0.3}, {'CM',3,0.1}}; %, {'Twist',2,0.2}};                
% QoI = 'Sectional_normal_force'; % Force at different radial stations
% %aero_module_outputfile = 'B1n_BEM.txt'; % Aero-Module data to be calibrated
% 
% % The parameter should be defined in the following format {name,index,rel_perturbation} where
% % rel_pertubation defines the amount of relative perturbation for B-spline
% % curves, and index the index associated with the B-spline control point.
% % For example, {'Twist',2,0.2} defines the uncertainty in the
% % second control point of Twist curve and 0.2 means a relative uncertainty of plus minus 10% for this  control point. 
% 
% % This parameter may not be required for other scalar random
% % variables; in that case we simply have {name,'',''} 
% 
% % uncertain_params = {{'Twist',2,0.2},{'Twist',3,0.2},{'Twist',4,0.2},{'Twist',5,0.2},{'Twist',6,0.2},{'Twist',7,0.2},...
% %                    {'Chord',2,0.2},{'Chord',4,0.2},{'Chord',6,0.2},{'Chord',8,0.2}, ...
% %                    {'Thickness',2,0.2},{'Thickness',3,0.2},{'Thickness',4,0.2},{'Thickness',5,0.2},...
% %                    {'YAW','',''},{'WINDSPEED','',''},{'RPM','',''},{'PITCHANGLE','',''}, ...
% %                    {'CL',1, 0.2}, {'CL',2, 0.2}, {'CL',3, 0.3},{'CL',4,0.3}, ...
% %                    {'CD',1, 0.2}, {'CD',2, 0.2}, {'CD',3, 0.2},{'CD',4,0.2}, ...
% %                    {'CM',1, 0.2}, {'CM',2, 0.2}, {'CM',3, 0.2},{'CM',4,0.2}
% %                    {'DYNSTALLTYPE','',''}, {'CORR3DTYPE','',''},...
% %                    {'BL_A1','',''},{'BL_A2','',''},{'BL_b1','',''},{'BL_b2','',''},...
% %                    {'BL_Ka','',''},{'BL_Tp','',''},{'BL_Tf','',''},{'BL_Tv','',''},{'BL_Tvl','',''},{'BL_Acd','',''}};
% 
% %% Variables of input file extracted from reference test case from DANAERO turbine NM80
% Inp.AEROMODEL = 1; % 1:BEM 2: AWSM
% Inp.TURBINETYPE = 1; %1:HAWT 2: VAWT
% 
% % radial positions used for twist, chord, thickness
% Inp.zB = [0.0250    0.0900    0.1650    0.2400    0.4650    0.6900    0.8150    0.9160    1.0150    1.1400    1.2650    1.3650    1.4650    1.5900    1.8150    1.9550    1.9830    2.0120    2.0400];
% 
% % Baseline chord
% Inp.ref_chord = [0.0900    0.0900    0.1650    0.2400    0.2070    0.1780    0.1660    0.1580    0.1500    0.1420    0.1340    0.1290    0.1230    0.1160    0.1020    0.0920    0.0820    0.0560    0.0110]; 
%    
% % Thickness by chord ratio 
% Inp.t_by_c = [1.0000    1.0000    0.6250    0.2500    0.2500    0.2500    0.2500    0.2300    0.2100    0.2100    0.2100    0.1960    0.1800    0.1800    0.1800    0.1800    0.1800    0.1800    0.1800]; 
%       
% % Baseline twist
% Inp.ref_twist = [0   0    8.2000   16.4000   12.1000    8.3000    7.1000    6.1000    5.5000    4.8000    4.0000    3.7000    3.2000    2.6000    1.5000    0.7000    0.4690    0.2310         0]; 
%      
% Inp.C14 = [-2.3000   -2.3000   -1.1000         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0         0];
%    
% Inp.xB = [     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];
% 
% Inp.yB = [     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0     0];
% 
% Inp.vectorLength = length(Inp.zB);
% 
% Inp.BLADELENGTH = 2.04;
% Inp.BLADEROOT = 0.21;
% 
% Inp.HUBHEIGHT = 5.49;
% Inp.PITCHANGLE = -2.3;
% 
% Inp.RPM = 425.1;
% Inp.TBEGIN = 0.0;
% Inp.TEND = 2.1206;
% Inp.TILTANGLE = 0.0;
% Inp.TIMESTEP = 0.00393;
% %[TIMESTEP,TEND] = wakepoints(RPM); % Routine to compute the TIMESTEP and  TEND using RPM
% Inp.XNAC2HUB = -2.13;
% Inp.YAWANGLE = 0.0;
% Inp.ZNAC2HUB = 0.37;
% 
% Inp.WINDSPEED = 15.01;
% Inp.DYNSTALLTYPE = 1; % 0: No DS 1:Snel1 2: Snel2 3:B-Leishmann 4: Onera
% Inp.CORR3DTYPE = 1;
% 
% Inp.AIRDENSITY  = 1.2049;
% Inp.TOWERBASERADIUS = 0.254;
% Inp.TOWERTOPRADIUS = 0.254;
% 
% Inp.NROFBEMELEMENTS = 14;
% 
% % BL model parameters
% Inp.BL_A1 = 0.30;
% Inp.BL_A2 = 0.70;
% Inp.BL_b1 = 0.14;
% Inp.BL_b2 = 0.53;
% Inp.BL_Ka = 0.75;
% Inp.BL_Tp = 1.50;
% Inp.BL_Tf = 5.00;
% Inp.BL_Tv = 6.00;
% Inp.BL_Tvl = 5.00;
% Inp.BL_Acd = 0.13;
% 
% % Yaw model parameters
% Inp.AM11 =	0.445;
% Inp.AM12 = -1.78;
% Inp.AM13 = 1.63;
% Inp.AM14 = -0.0543;
% Inp.AM15 = 0.367;
% Inp.AM21 = 0.0523;
% Inp.AM22 = -0.284;
% Inp.AM23 = 0.327;
% Inp.AM24 = -0.0134;
% Inp.AM25 = 0.144;
% Inp.PH11 = -51.2;
% Inp.PH12 = 1009.0;
% Inp.PH13 = -1383.0;
% Inp.PH14 = 387.0;
% Inp.PH15 = -260.0;
% Inp.PH21 = 296.0;
% Inp.PH22 = 60.9;
% 
% 
% %% Define properties of uncertain input in the UQLab format. 
% % We define this for all possible uncertain inputs and finally in the 
% % variable "uncertain_params" we specify which variables to be considered                                    
% % for the sensitivity analysis. NOTE that the ordering of variables is not important.
% 
% counter = 0; % To keep track of the number of uncertain variables. The value should be increase by one after adding every new random parameter. 
% 
% %% ===============TWIST====================
% % Using 12 control points for twist chosen heuristically. Variable "Input.Marginals.Index" stores the index of
% % control points. Here we assume same marginal for all control points.
% NRCP_TWIST = 12;
% for i=1:NRCP_TWIST
%     counter = counter+1;
%     Input.Marginals(counter).Name = 'Twist';
%     Input.Marginals(counter).Index = i;
%     Input.Marginals(counter).Type = 'Uniform'; 
%     Input.Marginals(counter).Parameters = [-0.5 0.5];
%     Input.Marginals(counter).Bounds = [-0.5 0.5];
% end
% 
% %% ===============CHORD====================
% % Using 9 control points for chord chosen heuristically. Variable "Input.Marginals.Index" stores the index of
% % control points.
% NRCP_CHORD = 9;
% for i=1:NRCP_CHORD
%     counter = counter+1;
%     Input.Marginals(counter).Name = 'Chord';
%     Input.Marginals(counter).Index = i;
%     Input.Marginals(counter).Type = 'Uniform'; 
%     Input.Marginals(counter).Parameters = [-0.5 0.5];
%     Input.Marginals(counter).Bounds = [-0.5 0.5];
% end
% 
% %% ===============THICKNESS====================
% % Using 9 control points for thickness chosen heuristically. Variable "Input.Marginals.Index" stores the index of
% % control points.
% NRCP_THICKNESS = 9;
% for i=1:NRCP_THICKNESS
%     counter = counter+1;
%     Input.Marginals(counter).Name = 'Thickness';
%     Input.Marginals(counter).Index = i;
%     Input.Marginals(counter).Type = 'Uniform'; 
%     Input.Marginals(counter).Parameters = [-0.5 0.5];
%     Input.Marginals(counter).Bounds = [-0.5 0.5];
% end
% 
% %% =======================YAW====================
% % Truncated Gaussian
% YAWANGLE   = Inp.YAWANGLE;
% YAW_Std = 2;  % Standard deviation
% YAW_LB = -10; % Lower bound of trucated Gaussian distribution
% YAW_UB = 10;  % Upper bound of trucated Gaussian distribution
% counter = counter+1;
% Input.Marginals(counter).Name = 'YAWANGLE';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [YAWANGLE, YAW_Std];
% Input.Marginals(counter).Bounds = [YAW_LB YAW_UB]; 
% 
% %% =======================PITCH==================
% % Truncated Gaussian
% PITCHANGLE   = Inp.PITCHANGLE;
% PITCHANGLE_Std = 1;  % Standard deviation
% PITCHANGLE_LB = -2; % Lower bound of trucated Gaussian distribution
% PITCHANGLE_UB = 2;  % Upper bound of trucated Gaussian distribution
% counter = counter+1;
% Input.Marginals(counter).Name = 'PITCHANGLE'; 
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PITCHANGLE, PITCHANGLE_Std];
% Input.Marginals(counter).Bounds = [PITCHANGLE_LB PITCHANGLE_UB];
% 
% %% =======================RPM====================
% % Truncated Gaussian
% RPM = Inp.RPM;
% RPM_Std = 1;  % Standard deviation
% RPM_LB = 10; % Lower bound of trucated Gaussian distribution
% RPM_UB = 14;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'RPM';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [RPM, RPM_Std];
% Input.Marginals(counter).Bounds = [RPM_LB RPM_UB]; 
% 
% %% =======================WINDSPEED==============
% WindSpeed_scale = 6.1;
% WindSpeed_shape = 50;
% counter = counter + 1;
% Input.Marginals(counter).Name = 'WINDSPEED';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Weibull'; 
% Input.Marginals(counter).Parameters = [WindSpeed_scale WindSpeed_shape]; % scale and shape parameter
% Input.Marginals(counter).Bounds = ''; % No bound needed for Weibull 
% 
% %% =======================DYNSTALLTYPE==============
% % Discrete variable with values 0,1,2,3,4 with 0: No DS 1:Snel1 2: Snel2
% % 3:B-Leishmann 4: Onera. We use uniform random variable from [0 5] sample models with equal probability
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'DYNSTALLTYPE';
% Input.Marginals(counter).Index = '';
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [0 20-10^-20]; % subtract a small bias as floor operator should result in 0,1,2,3,4 
% Input.Marginals(counter).Bounds = [0 20-10^-20];
% 
% 
% %% =======================CORR3DTYPE==============
% % Discrete variable with values 0,1  with 0: No correction 1: Snel. 
% % We use uniform random variable from [0 2] sample models with equal probability
% counter = counter+1;
% Input.Marginals(counter).Name = 'CORR3DTYPE';
% Input.Marginals(counter).Index = '';
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [0 8-10^-20]; % subtract a small bias as floor operator should result in 0,1
% Input.Marginals(counter).Bounds = [0 8-10^-20];
% 
% 
% %% ====================Polars=====================
% % Import the reference polar curves
% % using fullfile works both on Windows and Linux and Mac
% section1_file = fullfile('..','..','AEROmodule','NewMexico_calibrate','reference','risoA121.dat');
% % section2_file = fullfile('..','..','AEROmodule','NewMexico_calibrate','reference','cylinder.dat');
% section3_file = fullfile('..','..','AEROmodule','NewMexico_calibrate','reference','du91w250.dat');
% section4_file = fullfile('..','..','AEROmodule','NewMexico_calibrate','reference','naca64418_clean_Re0.7M.dat');
% [alpha1, CL_section1, CD_section1, CM_section1] = importPolars(section1_file);
% % [alpha2, CL_section2, CD_section2, CM_section2] = importPolars(section2_file);
% [alpha3, CL_section3, CD_section3, CM_section3] = importPolars(section3_file);
% [alpha4, CL_section4, CD_section4, CM_section4] = importPolars(section4_file);
% % set indices where the polars are perturbed
% alpha1_pert = find(alpha1<=20 & alpha1>=0);
% % alpha2_pert = find(alpha2<=20 & alpha2>=0);
% alpha3_pert = find(alpha3<=20 & alpha3>=0);
% alpha4_pert = find(alpha4<=20 & alpha4>=0);
% 
% % Important: the format of POLARS cell defined below should not be changed 
% Inp.POLARS = {3,... % Number of polar files
%     {'risoA121.dat', 'du91w250.dat', 'naca64418_clean_Re0.7M.dat'}, ... % Name of polar files
%     {'Riso A1-21', 'DU 91-W250', 'NACA 64418'}, ... % Airfoil_Name
%     {0.21, 0.25 , 0.18}, ... % thickness by chord ratio
%     {1.6E+06, 5.0E+05, 7.0E+05}, ... % Reynolds numbers
%     {CL_section1, CD_section1, CM_section1, alpha1, alpha1_pert}, ... % Cl, Cd, Cm data for section03
% ...%     {CL_section2, CD_section2, CM_section2, alpha2, alpha2_pert}, ...
%     {CL_section3, CD_section3, CM_section3, alpha3, alpha3_pert}, ...
%     {CL_section4, CD_section4, CM_section4, alpha4, alpha4_pert}, ...
%     };
% 
% % Define PDF for CL
% counter = counter+1;
% Input.Marginals(counter).Name = 'CL';
% Input.Marginals(counter).Airfoil = 'risoA121.dat';
% Input.Marginals(counter).AlphaPert = alpha1_pert;
% Input.Marginals(counter).Index = 1; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CL';
% Input.Marginals(counter).Airfoil = 'du91w250.dat';
% Input.Marginals(counter).AlphaPert = alpha3_pert;
% Input.Marginals(counter).Index = 3; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CL';
% Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
% Input.Marginals(counter).AlphaPert = alpha4_pert;
% Input.Marginals(counter).Index = 4; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];
% % 
% % counter = counter+1;
% % Input.Marginals(counter).Name = 'CL';
% % Input.Marginals(counter).Index = 4; % unique ID corresponding to airfoil
% % Input.Marginals(counter).Type = 'Uniform'; 
% % Input.Marginals(counter).Parameters = [-1.5 1.5];
% % Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% % Define PDF for CD
% counter = counter+1;
% Input.Marginals(counter).Name = 'CD';
% Input.Marginals(counter).Index = 1; % Corresponds to section 3
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CD';
% Input.Marginals(counter).Index = 2; % Corresponds to section 5
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CD';
% Input.Marginals(counter).Index = 3; % Corresponds to section 8
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CD';
% Input.Marginals(counter).Index = 4; % Corresponds to section 10
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% % Define PDF for CM
% counter = counter+1;
% Input.Marginals(counter).Name = 'CM';
% Input.Marginals(counter).Index = 1; % Corresponds to section 3
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CM';
% Input.Marginals(counter).Index = 2; % Corresponds to section 5
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CM';
% Input.Marginals(counter).Airfoil = 'du91w250.dat';
% Input.Marginals(counter).AlphaPert = alpha3_pert;
% Input.Marginals(counter).Index = 3; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CM';
% Input.Marginals(counter).Index = 4; % Corresponds to section 10
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-0.5 0.5];
% Input.Marginals(counter).Bounds = [-0.5 0.5];
% 
% 
% %% Inputs for BL model parameters
% % =======================BL_A1====================
% % Truncated Gaussian
% BL_A1   = Inp.BL_A1;
% BL_A1_Std = 1;  % Standard deviation
% BL_A1_LB = BL_A1 - BL_A1*0.1; % Lower bound of trucated Gaussian distribution
% BL_A1_UB = BL_A1 + BL_A1*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_A1';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_A1, BL_A1_Std];
% Input.Marginals(counter).Bounds = [BL_A1_LB BL_A1_UB]; 
% 
% % =======================BL_A2====================
% % Truncated Gaussian
% BL_A2   = Inp.BL_A2;
% BL_A2_Std = 1;  % Standard deviation
% BL_A2_LB = BL_A2 - BL_A2*0.1; % Lower bound of trucated Gaussian distribution
% BL_A2_UB = BL_A2 + BL_A2*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_A2';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_A2, BL_A2_Std];
% Input.Marginals(counter).Bounds = [BL_A2_LB BL_A2_UB]; 
% 
% % =======================BL_b1====================
% % Truncated Gaussian
% BL_b1   = Inp.BL_b1;
% BL_b1_Std = 1;  % Standard deviation
% BL_b1_LB = BL_b1 - BL_b1*0.1; % Lower bound of trucated Gaussian distribution
% BL_b1_UB = BL_b1 + BL_b1*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_b1';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_b1, BL_b1_Std];
% Input.Marginals(counter).Bounds = [BL_b1_LB BL_b1_UB];
% 
% % =======================BL_b2====================
% % Truncated Gaussian
% BL_b2   = Inp.BL_b2;
% BL_b2_Std = 1;  % Standard deviation
% BL_b2_LB = BL_b2 - BL_b2*0.1; % Lower bound of trucated Gaussian distribution
% BL_b2_UB = BL_b2 + BL_b2*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_b2';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_b2, BL_b2_Std];
% Input.Marginals(counter).Bounds = [BL_b2_LB BL_b2_UB];
% 
% % =======================BL_Ka====================
% % Truncated Gaussian
% BL_Ka   = Inp.BL_Ka;
% BL_Ka_Std = 1;  % Standard deviation
% BL_Ka_LB = BL_Ka - BL_Ka*0.1; % Lower bound of trucated Gaussian distribution
% BL_Ka_UB = BL_Ka + BL_Ka*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_Ka';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_Ka, BL_Ka_Std];
% Input.Marginals(counter).Bounds = [BL_Ka_LB BL_Ka_UB];
% 
% % =======================BL_Tp====================
% % Truncated Gaussian
% BL_Tp   = Inp.BL_Tp;
% BL_Tp_Std = 1;  % Standard deviation
% BL_Tp_LB = BL_Tp - BL_Tp*0.1; % Lower bound of trucated Gaussian distribution
% BL_Tp_UB = BL_Tp + BL_Tp*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_Tp';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_Tp, BL_Tp_Std];
% Input.Marginals(counter).Bounds = [BL_Tp_LB BL_Tp_UB];
% 
%  
% % =======================BL_Tf====================
% % Truncated Gaussian
% BL_Tf   = Inp.BL_Tf;
% BL_Tf_Std = 1;  % Standard deviation
% BL_Tf_LB = BL_Tf - BL_Tf*0.1; % Lower bound of trucated Gaussian distribution
% BL_Tf_UB = BL_Tf + BL_Tf*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_Tf';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_Tf, BL_Tf_Std];
% Input.Marginals(counter).Bounds = [BL_Tf_LB BL_Tf_UB];
% 
% % =======================BL_Tv====================
% % Truncated Gaussian
% BL_Tv   = Inp.BL_Tv;
% BL_Tv_Std = 1;  % Standard deviation
% BL_Tv_LB = BL_Tv - BL_Tv*0.1; % Lower bound of trucated Gaussian distribution
% BL_Tv_UB = BL_Tv + BL_Tv*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_Tv';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_Tv, BL_Tv_Std];
% Input.Marginals(counter).Bounds = [BL_Tv_LB BL_Tv_UB];
% 
% % =======================BL_Tvl====================
% % Truncated Gaussian
% BL_Tvl   = Inp.BL_Tvl;
% BL_Tvl_Std = 1;  % Standard deviation
% BL_Tvl_LB = BL_Tvl - BL_Tvl*0.1; % Lower bound of trucated Gaussian distribution
% BL_Tvl_UB = BL_Tvl + BL_Tvl*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_Tv1';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_Tvl, BL_Tvl_Std];
% Input.Marginals(counter).Bounds = [BL_Tvl_LB BL_Tvl_UB];
% 
% % =======================BL_Acd====================
% % Truncated Gaussian
% BL_Acd   = Inp.BL_Acd;
% BL_Acd_Std = 1;  % Standard deviation
% BL_Acd_LB = BL_Acd - BL_Acd*0.1; % Lower bound of trucated Gaussian distribution
% BL_Acd_UB = BL_Acd + BL_Acd*0.1;  % Upper bound of trucated Gaussian distribution 
% counter = counter+1;
% Input.Marginals(counter).Name = 'BL_Acd';
% Input.Marginals(counter).Index = ''; % Empty for scalar
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [BL_Acd, BL_Acd_Std];
% Input.Marginals(counter).Bounds = [BL_Acd_LB BL_Acd_UB];



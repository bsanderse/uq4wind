function Input  = example_turbine_file()


% below is a list of variables that could be considered uncertain

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
YAW_Mean = 0;
YAW_Std = 2;  % Standard deviation
YAW_LB = -10; % Lower bound of trucated Gaussian distribution
YAW_UB = 10;  % Upper bound of trucated Gaussian distribution
counter = counter+1;
Input.Marginals(counter).Name = 'YAWANGLE';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [YAW_Mean, YAW_Std];
Input.Marginals(counter).Bounds = [YAW_LB YAW_UB]; 

%% =======================PITCH==================
% Truncated Gaussian
PITCHANGLE_Mean   = 0;
PITCHANGLE_Std = 1;  % Standard deviation
PITCHANGLE_LB = -2; % Lower bound of trucated Gaussian distribution
PITCHANGLE_UB = 2;  % Upper bound of trucated Gaussian distribution
counter = counter+1;
Input.Marginals(counter).Name = 'PITCHANGLE'; 
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [PITCHANGLE_Mean, PITCHANGLE_Std];
Input.Marginals(counter).Bounds = [PITCHANGLE_LB PITCHANGLE_UB];

%% =======================RPM====================
% Truncated Gaussian
RPM_Mean = 0;
RPM_Std = 1;  % Standard deviation
RPM_LB = 10; % Lower bound of trucated Gaussian distribution
RPM_UB = 14;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'RPM';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [RPM_Mean, RPM_Std];
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
counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section03_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section05_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 2; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section08_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];% 

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section10_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];


counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section03_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section05_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 2; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section08_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];% 

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section10_ref.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-1.5 1.5];
Input.Marginals(counter).Bounds = [-1.5 1.5];


%% Inputs for BL model parameters
% =======================BL_A1====================
% Truncated Gaussian
BL_A1   = 0.3;
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
BL_A2   = 0.7;
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
BL_b1   = 0.14;
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
BL_b2   = 0.53;
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
BL_Ka   = 0.75;
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
BL_Tp   = 1.5;
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
BL_Tf   = 5.0;
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
BL_Tv   = 6.0;
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
BL_Tvl   = 5.0;
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
BL_Acd   = 0.13;
BL_Acd_Std = 1;  % Standard deviation
BL_Acd_LB = BL_Acd - BL_Acd*0.1; % Lower bound of trucated Gaussian distribution
BL_Acd_UB = BL_Acd + BL_Acd*0.1;  % Upper bound of trucated Gaussian distribution 
counter = counter+1;
Input.Marginals(counter).Name = 'BL_Acd';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [BL_Acd, BL_Acd_Std];
Input.Marginals(counter).Bounds = [BL_Acd_LB BL_Acd_UB];




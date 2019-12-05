function [AEROMODEL,TURBINETYPE,zB, ref_chord, t_by_c,ref_twist, C14, xB, yB, vectorLength, ...
          BLADELENGTH, BLADEROOT, HUBHEIGHT, TILTANGLE, PITCHANGLE, XNAC2HUB, ...
          RPM, TBEGIN, TEND, TIMESTEP, YAWANGLE, NROFBEMELEMENTS, ZNAC2HUB, Input, uncertain_params, QoI, WINDSPEED, POLARS]  = AVATAR()
%% Variables of input file extracted from reference test case from AVATAR turbine data
AEROMODEL = 1;
TURBINETYPE = 1;

zB = [0 0.6161 2.4491 5.4540 9.5568 14.6564 20.6272 27.3223 34.5768 42.2120 50.0400 ...
    57.8680 65.5032 72.7577 79.4528 85.4236 90.5232 94.6260 97.6309 99.4639 100.0800];
  
ref_chord = [5.3800 5.3800 5.3800 5.4216 5.6284 5.9043 6.1694 6.0687 5.7501 5.2967 ...
             4.8396 4.4002 3.9139 3.5039 3.1367 2.8092 2.5046 2.2438 1.3385 0.5601 0.2308]; % Baseline chord, not needed in the simulation...
     
t_by_c = [1.0000 0.9964 0.9896 0.9803 0.8791 0.6850 0.5066 0.4098 0.3597 0.3179 0.2787 ...
          0.2500 0.2462 0.2455 0.2450 0.2437 0.2426 0.2419 0.2406 0.2395 0.2392]; % Thickness by chord ratio should also be uncertain, currently considered deterministic
      
ref_twist = [17.2800 17.2800 17.2800 17.2800 16.6362 16.0567 13.9732 10.5030 8.6600 ...
             7.4956 6.6913 5.9479 5.1453 4.6420 4.2275 3.8807 3.6254 3.4425 3.2484 3.2132 3.2360]; % Baseline twist, not needed in the simulation...
     
C14 = [-25.0000 -24.9594 -24.9000 -24.4097 -21.6065 -17.7940 -12.9801 -10.2591 -10.0000 -10.0000 ...
       -10.0000 -10.0000 -10.0000 -10.0000 -10.0000 -10.0000 -10.0000 -10.0000 -10.0000 -10.0000 -10.0000];
   
xB = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

yB = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

vectorLength = length(zB);
BLADELENGTH = 100.08;
BLADEROOT = 2.8;
HUBHEIGHT = 132.7;
TILTANGLE = 0.0;
PITCHANGLE = 0.15;
XNAC2HUB = -7.1;
RPM = 9.0218;
TBEGIN =0;
TEND = 79.8067;
TIMESTEP = 0.036947542*2;
%[TIMESTEP,TEND] = wakepoints(RPM); % Routine to compute the TIMESTEP and  TEND using RPM
YAWANGLE = 0.0;
NROFBEMELEMENTS = 21;
ZNAC2HUB = 3.369;
WINDSPEED = 6.1;

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
YAW_Std = 2 ;  % Standard deviation
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
 

%% =======================RPM====================
% Truncated Gaussian
RPM_Std = 1;  % Standard deviation
RPM_LB = RPM-2; % Lower bound of trucated Gaussian distribution
RPM_UB = RPM+2;  % Upper bound of trucated Gaussian distribution
counter = counter+1;
Input.Marginals(counter).Name = 'RPM';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [RPM, RPM_Std];
Input.Marginals(counter).Bounds = [RPM_LB RPM_UB]; 
% % To visualize the above distribution, uncomment the following 
% hist(trandn(RPM_LB*ones(10000,1),RPM_UB*ones(10000,1)))

%% =======================WINDSPEED==============
WindSpeed_scale = WINDSPEED;
WindSpeed_shape = 50;
counter = counter + 1;
Input.Marginals(counter).Name = 'WINDSPEED';
Input.Marginals(counter).Index = ''; % Empty for scalar
Input.Marginals(counter).Type = 'Weibull'; 
Input.Marginals(counter).Parameters = [WindSpeed_scale WindSpeed_shape]; % scale and shape parameter
Input.Marginals(counter).Bounds = ''; % No bound needed for Weibull 

%% ====================Polars=====================
POLARS ={0};
% % To visualize above distribution, uncomment the following
% hist(wblrnd(WindSpeed_scale,WindSpeed_shape,[10000,1]),20)
 
%% Specify uncertain parameters to consider in the sensitivity analysis 
% The parameter should be defined in the following format {name,index,rel_perturbation} where
% rel_pertubation defines the amount of relative perturbation for B-spline curves

% uncertain_params = {{'Twist',2,0.2},{'Twist',3,0.2},{'Twist',4,0.2},{'Twist',5,0.2},{'Twist',6,0.2},{'Twist',7,0.2},...
%                 {'Chord',2,0.2},{'Chord',4,0.2},{'Chord',6,0.2},{'Chord',8,0.2}, ...
%                 {'Thickness',2,0.2},{'Thickness',3,0.2},{'Thickness',4,0.2},{'Thickness',5,0.2},...
%                 {'YAW','',''},{'WINDSPEED','',''},{'RPM','',''},{'PITCHANGLE','',''}};

uncertain_params = {{'YAW','',''},{'WINDSPEED','',''},{'RPM','',''},{'PITCHANGLE','',''}};

% Specify quantity of interest
QoI = 'Axial_Force'; % 'Axial_Force' or  'Power'




function [Parameters, UncertainInputs]  = getParameterAeroModule(turbineName)
% This routine returns the deterministic parameters stores in a matlab cell
% variable 'P'. The meaning of each element in 'P' is the described in the
% code. NOTE: The order of element in P should not be changed.

turbineData = str2func(turbineName);
% P ={}; 
%% Get turbine data and uncertain input specifications
% [AllMarginals, UncertainInputs, QoI, Parameters]  = turbineData();

UncertainInputs  = turbineData();


%% Other parameters related to AERO module software
% Top and the end of AeroPower.dat file
% Parameters.startRow = 2;
% Parameters.endRow = round(Parameters.TEND/Parameters.TIMESTEP)+ Parameters.startRow; % depends on the TEND and TIMESTEP
% Location where the AERO module is stored

%% Populate the elements of P. NOTE: The order of elements in P should not be changed!
% P = Inp;
% P{25}=Input;
% P{26}=uncertain_params;
% P{27}=QoI;
Parameters.turbineName = turbineName;

% 
% P{1} = Inp.AEROMODEL; 
% P{2} = Inp.TURBINETYPE;
% P{3} = Inp.zB; % the radial distance along the blade pitch axis from the blade root
% P{4} = Inp.ref_chord;
% P{5} = Inp.t_by_c;
% P{6} = Inp.ref_twist;
% P{7} = Inp.C14;
% P{8} = Inp.xB;
% P{9} = Inp.yB; 
% P{10} = Inp.vectorLength;
% P{11} = Inp.BLADELENGTH;
% P{12} = Inp.BLADEROOT;
% P{13} = Inp.HUBHEIGHT;
% P{14} = Inp.TILTANGLE;
% P{15} = Inp.PITCHANGLE;
% P{16} = Inp.XNAC2HUB;
% P{17} = Inp.RPM;
% P{18} = Inp.TEND;
% P{19} = Inp.TIMESTEP;
% P{20} = Inp.YAWANGLE;
% P{21} = Inp.NROFBEMELEMENTS; 
% P{22} = Inp.ZNAC2HUB;
% P{23} = Inp.startRow;
% P{24} = Inp.endRow;
% P{25} = Input;
% P{26} = uncertain_params;
% P{27} = QoI;
% P{28} = Inp.WINDSPEED;
% P{29} = turbineName;
% P{30} = Inp.TBEGIN;
% P{31} = Inp.POLARS;
% P{32} = Inp.DYNSTALLTYPE;
% P{33} = Inp.CORR3DTYPE;
% % BL model parameters
% P{34} = Inp.BL_A1;
% P{35} = Inp.BL_A2;
% P{36} = Inp.BL_b1;
% P{37} = Inp.BL_b2;
% P{38} = Inp.BL_Ka;
% P{39} = Inp.BL_Tp;
% P{40} = Inp.BL_Tf;
% P{41} = Inp.BL_Tv;
% P{42} = Inp.BL_Tvl;
% P{43} = Inp.BL_Acd;
% % parameters added when NewMexico case was introduced
% P{44} = Inp.AIRDENSITY;

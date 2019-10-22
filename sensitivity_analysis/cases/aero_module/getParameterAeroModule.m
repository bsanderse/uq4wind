function P = getParameterAeroModule(turbineName)
% This routine returns the deterministic parameters stores in a matlab cell
% variable 'P'. The meaning of each element in 'P' is the described in the
% code. NOTE: The order of element in P should not be changed.

turbineData = str2func(turbineName);
P ={}; 
%% Get turbine data and uncertain input specifications
[AEROMODEL,TURBINETYPE,zB, ref_chord, t_by_c,ref_twist, C14, xB, yB, vectorLength, ...
          BLADELENGTH, BLADEROOT, HUBHEIGHT, TILTANGLE, PITCHANGLE, XNAC2HUB, ...
          RPM, TBEGIN, TEND, TIMESTEP, YAWANGLE, NROFBEMELEMENTS, ZNAC2HUB, Input, uncertain_params, QoI, WINDSPEED]  = turbineData();

%% Other parameters related to AERO module software
% Top and the end of AeroPower.dat file
startRow = 2;
endRow = round(TEND/TIMESTEP)+ startRow; % depends on the TEND and TIMESTEP
% Location where the AERO module is stored
%% Populate the elements of P. The order of elements in P should not be changed!
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
P{25} = Input;
P{26} = uncertain_params;
P{27} = QoI;
P{28} = WINDSPEED;
P{29} = turbineName;
P{30} = TBEGIN;

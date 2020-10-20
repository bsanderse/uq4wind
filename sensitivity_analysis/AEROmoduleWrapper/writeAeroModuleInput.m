function writeAeroModuleInput(X,P)
% This routine add the random input to the input.txt file used for Aero
% module. Check getParameterAeroModule.m for order list of parameters
% stored in variable 'P'
%% ===========Get Twist samples===================
ndim = length(P{26});
TWIST_INDEX = [];
TWIST_PERTURB = [];
X_TWIST = [];
for i=1:ndim
    if(strcmp(P{26}{i}{1},'Twist'))
        TWIST_INDEX = [TWIST_INDEX P{26}{i}{2}];
        TWIST_PERTURB = [TWIST_PERTURB P{26}{i}{3}];
        X_TWIST = [X_TWIST X(i)];
    end
end

if(length(TWIST_PERTURB)>=1)
    % short hand:
    twist = computeTwist(1, TWIST_INDEX, X_TWIST, TWIST_PERTURB, 0); % computeTwist routine uses the specifications of NM80 turbine by default
    % full specification:
    %computeTwist(samples, index, randVec, pc, plotSamples,bladeLength,...
    %                                  interpolationLocations,referenceTwist,...
    %                                  t0,n,sampledLocations,sampledValues)
    %twist = computeTwist(1, TWIST_INDEX, X_TWIST, TWIST_PERTURB, 0, ...
    %     P{11},); % computeTwist routine uses the specifications of NM80 turbine by default
else
    twist = P{6};
end
%% ===========Get Chord samples===================
CHORD_INDEX = [];
CHORD_PERTURB = [];
X_CHORD = [];
for i=1:ndim
    if(strcmp(P{26}{i}{1},'Chord'))
        CHORD_INDEX = [CHORD_INDEX P{26}{i}{2}];
        CHORD_PERTURB = [CHORD_PERTURB P{26}{i}{3}];
        X_CHORD = [X_CHORD X(i)];
    end
end
if(length(CHORD_PERTURB)>=1)
    chord = computeChord(1,CHORD_INDEX, X_CHORD, CHORD_PERTURB, 0); % computeChord routine uses the specifications of NM80 turbine by default
else
    chord = P{4};
end
%% ===========Get Thickness samples===============
THICKNESS_INDEX = [];
THICKNESS_PERTURB = [];
X_THICKNESS = [];
for i=1:ndim
    if(strcmp(P{26}{i}{1},'Thickness'))
        THICKNESS_INDEX = [THICKNESS_INDEX P{26}{i}{2}];
        THICKNESS_PERTURB = [THICKNESS_PERTURB P{26}{i}{3}];
        X_THICKNESS = [X_THICKNESS X(i)];
    end
end
if(length(THICKNESS_INDEX)>=1)
    thickness = computeThickness(1, THICKNESS_INDEX, X_THICKNESS, THICKNESS_PERTURB, 0); % computeThickness routine uses the specifications of NM80 turbine by default
else
    thickness =P{5}.*P{4};
end

%% ===========Get YAW sample===============
X_YAW = P{20}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'YAW'))
        X_YAW = X(i);
    end
end

%% ===========DYNSTALLTYPE===============
X_DYNSTALLTYPE = P{32}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'DYNSTALLTYPE'))
%         X_DYNSTALLTYPE = floor(X(i)/4);
        X_DYNSTALLTYPE = X(i);
    end
end

%% =========== BL parameters ===============
X_BL_A1 = P{34}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_A1'))
        X_BL_A1 = X(i);
    end
end

X_BL_A2 = P{35}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_A2'))
        X_BL_A2 = X(i);
    end
end

X_BL_b1 = P{36}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_b1'))
        X_BL_b1 = X(i);
    end
end

X_BL_b2 = P{37}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_b2'))
        X_BL_b2 = X(i);
    end
end

X_BL_Ka = P{38}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_Ka'))
        X_BL_Ka = X(i);
    end
end

X_BL_Tp = P{39}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_Tp'))
        X_BL_Tp = X(i);
    end
end

X_BL_Tf = P{40}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_Tf'))
        X_BL_Tf = X(i);
    end
end

X_BL_Tv = P{41}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_Tv'))
        X_BL_Tv = X(i);
    end
end

X_BL_Tvl = P{42}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_Tvl'))
        X_BL_Tvl = X(i);
    end
end

X_BL_Acd = P{43}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'BL_Acd'))
        X_BL_Acd = X(i);
    end
end

%% ===========CORR3DTYPE===============
X_CORR3DTYPE = P{33}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'CORR3DTYPE'))
%         X_CORR3DTYPE = floor(X(i)/4);
        X_CORR3DTYPE = X(i);
    end
end

%% ===========Get WindSpeed sample===============
X_WindSpeed = P{28}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'WINDSPEED'))
        X_WindSpeed = X(i);
    end
end

%% ===========Get RPM sample===============
X_RPM = P{17}; % Assign nominal value
for i=1:ndim
    if(strcmp(P{26}{i}{1},'RPM'))
        X_RPM = X(i);
    end
end

%% ===========Get RPM sample===============
X_PITCHANGLE = P{15}; % Assign nominal value
for i = 1:ndim
    if(strcmp(P{26}{i}{1},'PITCHANGLE'))
        X_PITCHANGLE = X(i);
    end
end

%% ========== Get CL sample ==============
CL_INDEX = [];
CL_PERTURB = [];
X_CL = [];
n_polar = P{31}{1};
ind_aoa = P{31}{7+n_polar}; % indices that are perturbed
aoa  = P{31}{6};
% note the perturbed angle of attack is aoa(ind_aoa)

CL =cell(n_polar); % initialize, P{31} contains the polars
for i = 1:n_polar % Loop over all possible polar files, P{31}{1} contains number of polars
    CL{i} = P{31}{6+i}{1}; % CL of each polar, corresponding to a certain section
end

% note:
% P{name,index,rel_perturbation}
for i=1:ndim
    if(strcmp(P{26}{i}{1},'CL')) % check if CL is part of the random variable vector
        CL_INDEX = [CL_INDEX P{26}{i}{2}]; % index -> section / polar
        CL_PERTURB = [CL_PERTURB P{26}{i}{3}]; % perturbation 
        X_CL = [X_CL X(i)]; % value of the random variable CL
    end
end
% now that the value of Cl and perturbation are defined, 
% the Cl curves can be constructed
% i.e. the Cl-alpha curve at a certain section
d = ones(length(ind_aoa),1);
for i = 1:n_polar % Loop over all possible polar files
    for j = 1:length(CL_INDEX)
        if(CL_INDEX(j)==i)
            % note the syntax:
            % computeCurves(samples,index, randVec, pc, plotSamples,...
%                                       interpolationLocations,referenceCurve,NURBS_order,sampledIndices)
            % X_CL gives a sample between [-0.5,0.5]
            % this is multiplied by the value of CL_PERTURB
            % computeCurves constructs a NURBS in Cl-alpha space
            % ind_aoa consists of the indices that are perturbed
            plotCurve = 0;
            CL{i} = computeCurves(1, ind_aoa, X_CL(i)*d, CL_PERTURB(i)*d, plotCurve, ...
                aoa, P{31}{6+i}{1}, 3, 1:length(aoa));
        end
    end
end

%% ========== Get CD sample ==============
CD_INDEX = [];
CD_PERTURB = [];
X_CD = [];

CD =cell(n_polar); % initialize, P{31} contains the polars

for i = 1:n_polar % Loop over all possible polar files
    CD{i} = P{31}{6+i}{2}; % CD of each polar, corresponding to a certain section
end
for i=1:ndim
    if(strcmp(P{26}{i}{1},'CD'))
        CD_INDEX = [CD_INDEX P{26}{i}{2}];
        CD_PERTURB = [CD_PERTURB P{26}{i}{3}];
        X_CD = [X_CD X(i)];
    end
end
d = ones(length(ind_aoa),1);
for i = 1:n_polar % Loop over all possible polar files
    for j = 1:length(CD_INDEX)
        if(CD_INDEX(j)== i)
            plotCurve = 0;
            CD{i} = computeCurves(1, ind_aoa, X_CD(i)*d, CD_PERTURB(i)*d, plotCurve, ...
                aoa, P{31}{6+i}{2}, 3, 1:length(aoa));
%             CD{i} = computeCurves(1,P{31}{7+P{31}{1}}, X_CD(i)*ones(length(P{31}{7+P{31}{1}}),1), ...
%                 CD_PERTURB(i)*ones(length(P{31}{7+P{31}{1}}),1), 0, P{31}{6}, P{31}{6+i}{2},3,1:length(P{31}{6}));
        end
    end
end

%% ========== Get CM sample ==============
CM_INDEX = [];
CM_PERTURB = [];
X_CM = [];
CM = cell(n_polar);
for i = 1:n_polar % Loop over all possible polar files, P{31}{1} contains number of polars
    CM{i} = P{31}{6+i}{3}; % CM of each polar, corresponding to a certain section
end
for i=1:ndim
    if(strcmp(P{26}{i}{1},'CM'))
        CM_INDEX = [CM_INDEX P{26}{i}{2}];
        CM_PERTURB = [CM_PERTURB P{26}{i}{3}];
        X_CM = [X_CM X(i)];
    end
end
d = ones(length(ind_aoa),1);
for i = 1:n_polar % Loop over all possible polar files
    for j = 1:length(CM_INDEX)
        if(CM_INDEX(j)==i)
            plotCurve = 0;
            CM{i} = computeCurves(1, ind_aoa, X_CM(i)*d, CM_PERTURB(i)*d, plotCurve, ...
                aoa, P{31}{6+i}{3}, 3, 1:length(aoa));
%             CM{i} = computeCurves(1,P{31}{7+P{31}{1}}, X_CM(i)*ones(length(P{31}{7+P{31}{1}}),1), ...
%                 CM_PERTURB(i)*ones(length(P{31}{7+P{31}{1}}),1), 0, P{31}{6}, P{31}{6+i}{3}, 3, 1:length(P{31}{6}));
        end
    end
end

%% Write to the input.txt file for aeromodule
filename = fullfile(pwd,'AEROmodule',P{29},'current','input.txt');
%filename = [pwd,'\AEROmodule\',P{29},'\input.txt'];
fid = fopen(filename,'w');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! General ------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'AEROMODEL                        %d	! 1:BEM 2: AWSM\n', P{1}); 
fprintf(fid,'!TURBINETYPE                     %d	! 1:HAWT 2: VAWT\n', P{2});
fprintf(fid,'!INCLUDE                         specialist_input.txt\n');
fprintf(fid,'!INTERPOL                        2  	! 1:Linear 2: Spline\n');
fprintf(fid,'LOGFILENAME                      logfile.dat\n');
fprintf(fid,'DEBUGFILE                        1\n');
fprintf(fid,'!SUBITERFLAG                     0	! 0:dt>0 1:all dt 2:>=9+extr 3:1 and >=9\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Blade definition ---------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'AEROPROPS\n');
fprintf(fid,'!zB [m] chord [m] t/c [-] twist [deg]  C14 [%%c] xB [m] yB [m]\n');
% loop over all radial points
for i = 1:P{10}
    fprintf(fid,'%f    %f    %f    %f    %f    %f    %f \n', P{3}(i), chord(i), thickness(i)/chord(i), twist(i), P{7}(i), P{8}(i), P{9}(i));
end
fprintf(fid,'BLADELENGTH                     %f\n',P{11});
fprintf(fid,'BLADEROOT                       %f\n',P{12});
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Stand alone  -------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'CONEANGLE                       0.0\n');
fprintf(fid,'HUBHEIGHT                       %f\n',P{13});
fprintf(fid,'NROFBLADES                      3\n');     
fprintf(fid,'PITCHANGLE                      %f\n',X_PITCHANGLE);
fprintf(fid,'RPM                             %f\n',X_RPM);
fprintf(fid,'TBEGIN                          %f\n',P{30}); 
fprintf(fid,'TEND                            %6.15f		!~3D\n',P{18});
fprintf(fid,'TILTANGLE                       %f\n', P{14});
fprintf(fid,'TIMESTEP                        %6.15f	!10deg\n',P{19});
fprintf(fid,'XNAC2HUB                        %f\n', P{16});
fprintf(fid,'YAWANGLE                        %f\n', X_YAW);
fprintf(fid,'ZNAC2HUB                        %f\n', P{22});
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Airfoil data -------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');              
%fprintf(fid,'AOAEVAL                         1\n');
fprintf(fid,'COEFFFILENAME                   airfoils.dat\n');
fprintf(fid,'CORR3DTYPE            %d !0: No correction 1: Snel\n', X_CORR3DTYPE);
fprintf(fid,'DYNSTALLTYPE          %d!0: No DS 1:Snel1 2: Snel2 3:B-Leishmann 4: Onera\n', X_DYNSTALLTYPE);
fprintf(fid,'BL_A1	                    %f\n', X_BL_A1);
fprintf(fid,'BL_A2	                    %f\n', X_BL_A2);
fprintf(fid,'BL_B1	                    %f\n', X_BL_b1);
fprintf(fid,'BL_B2	                    %f\n', X_BL_b2);
fprintf(fid,'BL_KA	                    %f\n', X_BL_Ka);
fprintf(fid,'BL_TP	                    %f\n', X_BL_Tp);
fprintf(fid,'BL_TF	                    %f\n', X_BL_Tf);
fprintf(fid,'BL_TV	                    %f\n', X_BL_Tv);
fprintf(fid,'BL_TVL	                    %f\n', X_BL_Tvl);
fprintf(fid,'BL_ACD	                    %f\n', X_BL_Acd);
fprintf(fid,'FSMETHOD	                2\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Environment --------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'AIRDENSITY                      1.231\n');
fprintf(fid,'HORSHEAR                        0.0\n');
fprintf(fid,'SHEAREXP                        0.0\n');
fprintf(fid,'DYNVISC                         1.7879e-5\n');
%fprintf(fid,'HORSHEAR                        0.0\n');
%fprintf(fid,'SOS                             340.3\n'); 
%fprintf(fid,'RAMPTIME               		    0\n');
%fprintf(fid,'RAMPFACTOR                      1\n');
fprintf(fid,'TOWERBASERADIUS                 0.0001\n');
fprintf(fid,'TOWERTOPRADIUS                  0.0001\n');
fprintf(fid,'WINDFILENAME                    wind.dat\n');
fprintf(fid,'WINDTYPE                        1	! 1: Simple 2: SWIFT\n');
%fprintf(fid,'Z0                              0\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'!BEM  ----------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
% fprintf(fid,'!AEROROOT                        0.40\n');	
% fprintf(fid,'ATRANSITION                     0.38\n');                           
% fprintf(fid,'DEBUGFILE                       1\n');
% fprintf(fid,'DYNINFLOW                       1\n');
% fprintf(fid,'INDUCTIONDRAG                   0\n');
fprintf(fid,'NROFBEMELEMENTS                   %d\n', P{21});
% fprintf(fid,'PRANDTLROOT                     1\n');             
% fprintf(fid,'PRANDTLTIP                      1\n');                                         
% fprintf(fid,'!---------------------------------------------------------------------\n');
% fprintf(fid,'!AWSM ----------------------------------------------------------------\n');
% fprintf(fid,'!---------------------------------------------------------------------\n');
% fprintf(fid,'AERORESULTSFILE                 Aeroresultsout.dat\n');
% fprintf(fid,'GEOMRESULTSFILE                 geomresultsout.dat\n');
% fprintf(fid,'TIMESTARTAERO                   1\n');
% fprintf(fid,'TIMEENDAERO                     1440\n');
% fprintf(fid,'AEROOUTPUTINCREMENT             1\n');
% fprintf(fid,'TIMESTARTGEOM                   1\n');
% fprintf(fid,'TIMEENDGEOM                     1440\n');
% fprintf(fid,'GEOMOUTPUTINCREMENT             1\n');
% fprintf(fid,'TRACELEVEL                      2\n');
% fprintf(fid,'STDOUTSWITCH                    1   !0: no 1: yes\n');
% fprintf(fid,'LOGTRACESWITCH                  1   !0: no 1: yes\n');
% fprintf(fid,'CENTEREDCONTRPOINTS             0   !0: no 1: yes\n');
% fprintf(fid,'STREAMWISEWAKEPOINTS            1440\n');
% fprintf(fid,'FREESTRMWISEWAKEPOINTS          1440\n');
% fprintf(fid,'VORTEXCUTOFF                    2\n');
% fprintf(fid,'WAKECUTOFFRADIUS                0.2\n');
% fprintf(fid,'CIRCCONVCRIT                    0.0001\n');
% fprintf(fid,'LIFTCUTOFFRADIUS                0.01\n');    
% fprintf(fid,'NEWRAPHNEIGHBOURS               2\n');
% fprintf(fid,'STARTSMOOTHINGS                 3\n');
% fprintf(fid,'NUMITEQMAX                      10\n');
% fprintf(fid,'NUMITAEROCONV                   100\n');
% fprintf(fid,'EXTFIELDFLAG                    0   !0: off 1: on\n');
% fprintf(fid,'GROUNDFLAG                      0   !0: off 1: on\n');
% fprintf(fid,'PRSCRBWAKE               	    0\n');
% fprintf(fid,'CONVECFACTOR			        0\n');
fclose(fid);

%% Change windspeed file
filename = fullfile(pwd,'AEROmodule',P{29},'current','wind.dat');
% filename = [pwd,'\AEROmodule\',P{29},'\wind.dat'];
fid = fopen(filename,'w');
fprintf(fid,'!time [s]  u [m/s]  v [m/s]  w [m/s]\n');
fprintf(fid,'%f    %f    %f    %f\n', 0.0,  X_WindSpeed, 0.0, 0.0);
fprintf(fid,'%f    %f    %f    %f\n', 5.0,  X_WindSpeed, 0.0, 0.0);
fprintf(fid,'%f    %f    %f    %f\n', 10.0, X_WindSpeed, 0.0, 0.0);
fprintf(fid,'%f    %f    %f    %f\n', 20.0, X_WindSpeed, 0.0, 0.0);
fclose(fid);

%% Change the Polar file
for i = 1:P{31}{1} % Loop over the polar files
    filename = fullfile(pwd,'AEROmodule',P{29},'current',P{31}{2}{i});
%     filename = [pwd,'\AEROmodule\',P{29},'\',P{31}{2}{i}]; 
    fid = fopen(filename,'w');
    fprintf(fid,'! Aero module input file for airfoil data\n');
    fprintf(fid,'\n');
    fprintf(fid, 'Airfoil_Name  %s\n',P{31}{3}{i});
    fprintf(fid,'t/c            %f         ! thickness ratio w.r.t. chord\n', P{31}{4}{i});
    fprintf(fid,'\n');
    fprintf(fid,'format 1       !  1: alfa-cl-cd-cm	; 2: alfa-cl; alfa-cd; alfa-cm\n');
    fprintf(fid,'\n');
    fprintf(fid,'Reynolds_Nr %f\n',P{31}{5});
    % loop over angle of attack
    for j = 1:length(P{31}{6})
        fprintf(fid,'%f    %f    %f    %f\n', P{31}{6}(j), CL{i}(j), CD{i}(j), CM{i}(j));
    end
    fclose(fid);
end


%% Change specialist.txt file
filename = fullfile(pwd,'AEROmodule',P{29},'current','specialist_input.txt');
% filename = [pwd,'\AEROmodule\',P{29},'\specialist_input.txt'];
fid = fopen(filename,'w');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! General ------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'!TURBINETYPE                     1	! 1:HAWT 2: VAWT\n');
fprintf(fid,'!INTERPOL                        2  	! 1:Linear 2: Spline\n');
fprintf(fid,'!SUBITERFLAG                     1	! 0:dt>0 1:all dt 2:>=9+extr 3:1 and >=9\n');
fprintf(fid,'DEBUGFILE                       1\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Airfoil data ------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'!AOAEVAL                         1\n');
fprintf(fid,'!LINREG 	                 1\n');
fprintf(fid,'!MINAOALIN 	           	-3.000	!-5.00\n');
fprintf(fid,'!MAXAOALIN 	            	 5.000	!8.00\n');
fprintf(fid,'!---\n');
fprintf(fid,'!Beddoes Leishman Parameters -----------------------------------------\n');
fprintf(fid,'!---\n');
% fprintf(fid,'BL_A1	                    0.30\n');
% fprintf(fid,'BL_A2	                    0.70\n');
% fprintf(fid,'BL_b1	                    0.14\n');
% fprintf(fid,'BL_b2	                    0.53\n');
% fprintf(fid,'BL_Ka	                    0.75\n');
% fprintf(fid,'BL_Tp	                    1.50\n');
% fprintf(fid,'BL_Tf	                    5.00\n');
% fprintf(fid,'BL_Tv	                    6.00\n');
% fprintf(fid,'BL_Tvl	                    5.00\n');
% fprintf(fid,'BL_Acd	                    0.13\n');
fprintf(fid,'BL_A1	                    %f\n', X_BL_A1);
fprintf(fid,'BL_A2	                    %f\n', X_BL_A2);
fprintf(fid,'BL_B1	                    %f\n', X_BL_b1);
fprintf(fid,'BL_B2	                    %f\n', X_BL_b2);
fprintf(fid,'BL_KA	                    %f\n', X_BL_Ka);
fprintf(fid,'BL_TP	                    %f\n', X_BL_Tp);
fprintf(fid,'BL_TF	                    %f\n', X_BL_Tf);
fprintf(fid,'BL_TV	                    %f\n', X_BL_Tv);
fprintf(fid,'BL_TVL	                    %f\n', X_BL_Tvl);
fprintf(fid,'BL_ACD	                    %f\n', X_BL_Acd);
fprintf(fid,'FSMETHOD	                2\n');
fprintf(fid,'!---\n');
fprintf(fid,'!ONERA Parameters -----------------------------------------\n');
fprintf(fid,'!---\n');
fprintf(fid,'!ON_LAML	                    0.17\n');
fprintf(fid,'!ON_SIGL	                    6.28\n');
fprintf(fid,'!ON_R0  	                    0.20\n');
fprintf(fid,'!ON_R2  	                    0.20\n');
fprintf(fid,'!ON_A0  	                    0.30\n');
fprintf(fid,'!ON_A2  	                    0.20\n');
fprintf(fid,'!ON_E2  	                    0.53\n');
fprintf(fid,'!ON_SL  	                    3.14\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Environment ------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'!SOS                             340.3 \n');
fprintf(fid,'!RAMPTIME			0.0\n');
fprintf(fid,'!RAMPFACTOR			1.0\n');
fprintf(fid,'!WINDINTERPL			2	! 1:Linear 2: Cubic (Turbsim/Mann only)\n');
fprintf(fid,'!Z0                              0.0\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! BEM ------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'!AEROROOT                        0\n');
fprintf(fid,'!ATRANSITION                     0.38  \n');
fprintf(fid,'!COUNTERMIN			50\n');
fprintf(fid,'!DYNINFLOW                       1\n');
fprintf(fid,'!INDUCTIONDRAG                   0\n');
fprintf(fid,'!PRANDTLROOT                     1 \n');
fprintf(fid,'!PRANDTLTIP                      1\n');
fprintf(fid,'!YAWMODEL			1\n');
fprintf(fid,'!SOLVEBEMMETHOD			1\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! AWSM ------------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'!AERORESULTSFILE                 aeroresultsout.dat\n');
fprintf(fid,'!GEOMRESULTSFILE                 geomresultsout.dat\n');
fprintf(fid,'!TIMESTARTAERO                   1\n');
fprintf(fid,'!TIMEENDAERO                     540\n');
fprintf(fid,'!AEROOUTPUTINCREMENT             1\n');
fprintf(fid,'!TIMESTARTGEOM                   540\n');
fprintf(fid,'!TIMEENDGEOM                     540\n');
fprintf(fid,'!GEOMOUTPUTINCREMENT             1\n');
fprintf(fid,'!TRACELEVEL                      2\n');
fprintf(fid,'!STDOUTSWITCH                    0   !0: no 1: yes\n');
fprintf(fid,'!LOGTRACESWITCH                  0   !0: no 1: yes\n');
fprintf(fid,'!CENTEREDCONTRPOINTS             0   !0: no 1: yes\n');
fprintf(fid,'!VORTEXCUTOFF                    2\n');
fprintf(fid,'!WAKECUTOFFRADIUS                0.1\n');
fprintf(fid,'!CIRCCONVCRIT                    0.0001\n');
fprintf(fid,'!LIFTCUTOFFRADIUS                0.001\n');
fprintf(fid,'!NEWRAPHNEIGHBOURS               2\n');
fprintf(fid,'!STARTSMOOTHINGS                 3\n');
fprintf(fid,'!NUMITEQMAX                      10\n');
fprintf(fid,'!NUMITAEROCONV                   100\n');
fprintf(fid,'!EXTFIELDFLAG                    0   !0: off 1: on\n');
fprintf(fid,'!GROUNDFLAG                      0   !0: off 1: on\n');
fprintf(fid,'!GROUNDLEVEL	                0.0\n');
fprintf(fid,'!PARALLEL	                1   !0: off 1: on\n');
fprintf(fid,'!PRSCRBWAKE			0   !0: off 1: on\n');
fprintf(fid,'!CONVECFACTOR			0.0\n');
fprintf(fid,'!WAKEREDUCTIONSTART		360\n');
fprintf(fid,'!WAKEREDUCTIONSKIP		4\n');
fclose(fid);

%% Uncomment to plot the random samples of chord, twist, CL, CD
% figure(1)
% hold on
% plot(P{3}, chord, 'color','k','linestyle','--','linewidth',0.5,'HandleVisibility','off')
% 
% figure(2)
% hold on 
% plot(P{3}, twist, 'color','k','linestyle','--','linewidth',0.5,'HandleVisibility','off')
% 
% figure(3)
% hold on
% plot(P{31}{6}, CL{2},  'color','k','linestyle','--','linewidth',0.5,'HandleVisibility','off')
% 
% figure(4)
% hold on 
% plot( P{31}{6}, CD{2},'color','k','linestyle','--','linewidth',0.5,'HandleVisibility','off')

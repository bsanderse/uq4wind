function writeAeroModuleInput(X,P)
% This routine add the random input to the input.txt file used for Aero module. 
%% ===========Get Twist samples===================
ndim = length(P{26});
TWIST_INDEX = [];
X_TWIST = [];
for i=1:ndim
    if(strcmp(P{26}{i}{1},'Twist'))
        TWIST_INDEX = [TWIST_INDEX P{26}{i}{2}];
        X_TWIST = [X_TWIST X(i)];
    end
end
twist = computeTwist(1, X_TWIST, 0.2*ones(X_TWIST), 0); % computeTwist routine uses the specifications of NM80 turbine by default

%% ===========Get Chord samples===================
CHORD_INDEX = [];
X_CHORD = [];
for i=1:ndim
    if(strcmp(P{26}{i}{1},'Chord'))
        CHORD_INDEX = [CHORD_INDEX P{26}{i}{2}];
        X_CHORD = [X_CHORD X(i)];
    end
end
chord = computeChord(1, X_CHORD, 0.2*ones(X_CHORD), 0); % computeTwist routine uses the specifications of NM80 turbine by default

%% ===========Get Thickness samples===============
THICKNESS_INDEX = [];
X_THICKNESS = [];
for i=1:ndim
    if(strcmp(P{26}{i}{1},'Thicknes'))
        THICKNESS_INDEX = [THICKNESS_INDEX P{26}{i}{2}];
        X_THICKNESS = [X_THICKNESS X(i)];
    end
end
thickness = computeThickness(1, X_THICKNESS, 0.2*ones(X_THICKNESS), 0); % computeTwist routine uses the specifications of NM80 turbine by default

                          
filename = [folder,'input.txt'];
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
for i = 1:vectorLength
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
fprintf(fid,'PITCHANGLE                      %f\n',P{15});
fprintf(fid,'RPM                             %f\n',P{17});
fprintf(fid,'TBEGIN                          0.0\n');
fprintf(fid,'TEND                            %6.15f		!~3D\n',P{18});
fprintf(fid,'TILTANGLE                       %f\n', P{14});
fprintf(fid,'TIMESTEP                        %6.15f	!10deg\n',P{19});
fprintf(fid,'XNAC2HUB                        %f\n', P{16});
fprintf(fid,'YAWANGLE                        %f\n', P{20});
fprintf(fid,'ZNAC2HUB                        %f\n', P{22});
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Airfoil data -------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');              
%fprintf(fid,'AOAEVAL                         1\n');
fprintf(fid,'COEFFFILENAME                   airfoils.dat\n');
fprintf(fid,'CORR3DTYPE            0 !0: No correction 1: Snel\n');
fprintf(fid,'DYNSTALLTYPE          1 !0: No DS 1:Snel1 2: Snel2 3:B-Leishmann 4: Onera\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'! Environment --------------------------------------------------------\n');
fprintf(fid,'!---------------------------------------------------------------------\n');
fprintf(fid,'AIRDENSITY                      1.231\n');
fprintf(fid,'HORSHEAR                        0.0\n');
fprintf(fid,'SHEAREXP                        0.0\n');
%fprintf(fid,'DYNVISC                         1.7879e-5\n');
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
end
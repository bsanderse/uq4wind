function statis(daqwin)

name = 'Wind speed at 17m (Wsp_17m) [m/s]';
signalnumber = 118;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind speed at 28.5m (Wsp_28.5m) [m/s]';
signalnumber = 117;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind speed at 41m (Wsp_41m) [m/s]';
signalnumber = 116;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind speed at 57m (Wsp_57m) [m/s]';
signalnumber = 115;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind speed at 77m (Wsp_77m) [m/s]';
signalnumber = 114;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind speed at 93m (Wsp_90m) [m/s]';
signalnumber = 113;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

%

name = 'Wind direction at 17m (Wind_dir_17m) [deg]';
signalnumber = 132;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind direction at 57m (Wind_dir_57m) [deg]';
signalnumber = 131;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Wind direction at 93m (Wind_dir_90m) [deg]';
signalnumber = 130;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

%

name = 'Yaw position (Yaw_Pos) [deg]';
signalnumber = 70;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Yaw error (-) [deg]';
signalnumber = 175;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

%

name = 'Fx03 (Fx03) [N/m]';
signalnumber = 138;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Fx05 (Fx05) [N/m]';
signalnumber = 139;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Fx08 (Fx08) [N/m]';
signalnumber = 140;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Fx10 (Fx10) [N/m]';
signalnumber = 141;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

%

name = 'Fy03 (Fy03) [N/m]';
signalnumber = 142;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Fy05 (Fy05) [N/m]';
signalnumber = 143;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Fy08 (Fy08) [N/m]';
signalnumber = 144;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Fy10 (Fy10) [N/m]';
signalnumber = 145;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

%

name = 'Power (Pow) [kW]';
signalnumber = 83;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Rotation speed (Rot_spd) [rpm]';
signalnumber = 84;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Pitch (Pitch_v1) [deg]';
signalnumber = 86;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

name = 'Density (rho) [kg/m^3]';
signalnumber = 137;
[meanout,stdout,minout,maxout] = figstat(daqwin,name,signalnumber);

end
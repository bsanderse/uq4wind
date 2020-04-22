%%
clc
close all
clearvars
format compact

%% Independent path
addpath([pwd,'/AEROmoduleWrapper/']);

%% Read AeroModule data
%run('AEROmodule/NM80_calibrate/output/read_interpolate.m');
filename_aero = ('AEROmodule/NM80_calibrate/output\output.dat'); % location of the AeroPower.dat output file
output_a = importfile(filename_aero, 2);
axial_aero = output_a.aero_data;
rad = output_a.r_i;

%%
%run('../Experimental/WINDTRUE/main_exp.m');
filename_exp = ('../Experimental/WINDTRUE\output_e.dat');
output_e = importfile1(filename_exp, 2);
axial_exp = output_e.exp_data;

%%
plot(rad,axial_aero,'g--o')
hold on
plot(rad,axial_exp,'r--x')
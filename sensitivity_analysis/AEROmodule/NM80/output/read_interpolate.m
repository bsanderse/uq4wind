
%%
clc
close all
clearvars
format compact

%% Read data
% The AeroModule output files are stored in the following format:
% for e.g: B1n_BEM.txt
% where, B1 is the Blade number 1
% n is the normal force (Variable of interest)
% BEM indicates the calculations performed using Blade Element Momentum
% theory used in the AeroModule code.

% For choosing any variable of interest, simply replace (Ctrl+F) the blade
% number and the variable of interest in the script (Replace all)
% for e.g. B1n => B3vInd (Replace all)

B1n = readtable("B1n_BEM.txt",'HeaderLines', 4,"ReadVariableNames",true,...
    "PreserveVariableNames",true); % Reads the specific .txt file

%% Mean
% The mean values at the available (a) radial points are stored in a1..a26
% The values corresponds to the radial locations stored in r_a

a1 = mean(B1n.("0.1833"));
a2 = mean(B1n.("1.7549"));
a3 = mean(B1n.("3.3294"));
a4 = mean(B1n.("4.9066"));
a5 = mean(B1n.("6.4848"));
a6 = mean(B1n.("8.0632"));
a7 = mean(B1n.("9.6423"));
a8 = mean(B1n.("11.2252"));
a9 = mean(B1n.("12.8089"));
a10 = mean(B1n.("14.3906"));
a11 = mean(B1n.("15.9682"));
a12 = mean(B1n.("17.5466"));
a13 = mean(B1n.("19.1253"));
a14 = mean(B1n.("20.7040"));
a15 = mean(B1n.("22.2825"));
a16 = mean(B1n.("23.8605"));
a17 = mean(B1n.("25.4378"));
a18 = mean(B1n.("27.0144"));
a19 = mean(B1n.("28.5901"));
a20 = mean(B1n.("30.1649"));
a21 = mean(B1n.("31.7388"));
a22 = mean(B1n.("33.3118"));
a23 = mean(B1n.("34.8842"));
a24 = mean(B1n.("36.4562"));
a25 = mean(B1n.("37.6335"));
a26 = mean(B1n.("38.4193"));

B1n_a = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,...
    a20,a21,a22,a23,a24,a25,a26];
r_a = [0.1833,1.7549,3.3294,4.9066,6.4848,8.0632,9.6423,11.2252,12.8089,14.3906...
    15.9682,17.5466,19.1253,20.7040,22.2825,23.8605,25.4378,27.0144,28.5901...
    30.1649,31.7388,33.3118,34.8842,36.4562,37.6335,38.4193];



%% Interpolation

% The radial locations at which the measurements (experiments) are
% available to be obtained using spline interpolation.

r_i = [13,19,30,37]; % Measurement locations
B1n_i = spline(r_a,B1n_a,r_i); % Interpolated variable values

%% Plotting

plot(r_a,B1n_a,"ko")
hold on
plot(r_i,B1n_i,"rx",'markersize',10)
xlabel('r [m]')
ylabel('B1n')
legend('available','interpolated')
title('NM80 turbine')

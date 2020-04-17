
%%
clc
close all
clearvars
format compact

%% Read data
filename = 'B1n_BEM.txt'; % Specify the .txt file (Variable of interest)
D = readtable(filename,'HeaderLines', 4,"ReadVariableNames",true,...
    "PreserveVariableNames",true); % Reads the variable data from the 
                                   % specified .txt file
r_a = str2double(D.Properties.VariableNames(3:end)); % Radial stations
%% Mean
D_a = mean(D{:,3:end},1); % Mean values at different radial stations

%% Interpolation
r_i = [13,19,30,37]; % Measurement radial stations
D_i = spline(r_a,D_a,r_i); % Interpolated data

%% Write interpolated data
T = table(r_i',D_i');
writetable(T,'output.txt','Delimiter','\t','WriteRowNames',true);
type output.txt





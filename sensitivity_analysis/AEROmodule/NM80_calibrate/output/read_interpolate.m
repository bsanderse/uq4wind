function data = read_interpolate(filename)
%% Read output data generated by AeroModule
D = readtable(filename,'HeaderLines', 4,"ReadVariableNames",true,...
    "PreserveVariableNames",true); % Reads the variable data from the 
                                   % specified .txt file
r_a = str2double(D.Properties.VariableNames(3:end)); % Radial stations
%% Calculate mean
D_a = mean(D{:,3:end},1); % Mean values at different radial stations
%% Interpolation
% These values are made available from Danaero experiments
r_i = [13,19,30,37]; % Measurement radial stations
data = spline(r_a,D_a,r_i); % Interpolated data using spline



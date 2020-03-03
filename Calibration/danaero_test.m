clc
close all
clearvars

%% DANAERO experiment
N_data = 20;  % Number of available dataset
F_max = 0.01; % Maximum azimuthally averaged spanwise force (N)
% Normal force (experimental data) in Newtons
F_n = [0.27; 0.28; 0.29; 0.3; 0.31 ; 0.32; 0.33;  0.53; 0.54; 0.55; ...
        0.56; 0.57; 0.80; 0.81; 0.82; 0.83; 0.84; 0.86; 0.87; 0.88]*F_max;
%From BEM code
% Beta angle distribution
beta_r = [5; 10; 15; 17]*(pi/180); % in radians
beta = repelem(beta_r,[7 5 5 3]); % assuming constant across distribution
% Chord distribution
c_r = [0.9; 0.6; 0.4; 0.3];
c = repelem(c_r,[7 5 5 3]); % assuming constant across distribution
% Relative velocity distribution
v_r = [22; 35; 50; 56];
v = repelem(v_r,[7 5 5 3]); % assuming constant across distribution

% DANAERO plot (modified)
figure
plot(beta,F_n,'.')
xlabel('\beta [radians]')
ylabel('F_{n} [-]')
title('Modified experimental data')

%% Frequentist approach
format compact
format long
p = 2; % Number of parameters to calibrate
A = zeros(N_data,p); % size of design matrix, normally N_data>p
q_test = eye(p,p);
for i=1:p
    A(:,i) = dan_model(q_test(i,:),beta);
end

q_freq = regress(F_n,A)
F_n_OLS = dan_model(q_freq,beta);
cl = (q_freq(1)./0.5*1.225.*v_r.*c_r)
cd = (q_freq(2)./0.5*1.225.*v_r.*c_r)
figure
plot(beta, F_n, 'kx')
hold on
plot(beta,F_n_OLS, 'g--x')
xlabel('\beta')
ylabel('F_{n}')
title('Experiment vs. Frequentist')
legend('measurement data','OLS');


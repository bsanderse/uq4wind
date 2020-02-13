
% Disease prevelance example
clear; close all; clc;

%% Prior parameters 
a = 1;
b = 4;

%% Data
N = 1000;
X = 40;

%% Model
theta = linspace(0,1,100);
Y_prior = betapdf(theta,a,b);
Y_likelihood = nchoosek(N,X)*(theta.^X).*(1-theta).^(N-X);
Y_posterior = betapdf(theta,a+X,b+N-X);

%% Plots
h = figure(1);
subplot(3,1,1),
plot(theta,Y_prior,'LineWidth',3)
title('Prior','FontSize',20)
set(gca,'FontSize',20)
ylabel('pdf')
subplot(3,1,2),
plot(theta,Y_likelihood,'m','LineWidth',3)
set(gca,'FontSize',20)
title('Likelihood','FontSize',20)
ylabel('likelihood')
subplot(3,1,3),
plot(theta,Y_posterior,'r','LineWidth',3)
title('Posterior','FontSize',20)
ylabel('pdf')
set(h,'Position',[1000 150 900 900])
set(gca,'FontSize',20)
xlabel('Theta')
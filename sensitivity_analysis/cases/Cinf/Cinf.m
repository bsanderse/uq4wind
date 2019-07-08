function Y = Cinf(X)
% Cinf is an infinitely differentiable function
% see Patterson,
% ﻿On the construction of a practical Ermakov-Zolotukhin multiple integrator
% integral 5

% computing the response value
Y(:,1) = cos(0.6*pi + 0.5*sum(X,2)); 

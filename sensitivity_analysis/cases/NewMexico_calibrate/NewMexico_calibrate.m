function Input  = NewMexico_calibrate()

% initialize counter of marginals
counter = 0;

% define marginals of yaw model
counter = counter+1;
AM11     = 0.446;
AM11_Std = 0.1;  % Standard deviation
AM11_LB = 0; % Lower bound of truncated Gaussian distribution
AM11_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM11';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM11, AM11_Std];
Input.Marginals(counter).Bounds = [AM11_LB AM11_UB]; 

% define marginals 
% counter = counter+1;
% Input.Marginals(counter).Name = 'CL';
% Input.Marginals(counter).Airfoil = 'risoA121.dat';
% Input.Marginals(counter).AlphaPert = [-10 30];
% Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];

% counter = counter+1;
% Input.Marginals(counter).Name = 'CL';
% Input.Marginals(counter).Airfoil = 'du91w250.dat';
% Input.Marginals(counter).AlphaPert = [-10 30];
% Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];

% counter = counter+1;
% Input.Marginals(counter).Name = 'CL';
% Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
% Input.Marginals(counter).AlphaPert = [-10 30];
% Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];% 

% counter = counter+1;
% Input.Marginals(counter).Name = 'CD';
% Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
% Input.Marginals(counter).AlphaPert = [0 20];
% Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
% Input.Marginals(counter).Type = 'Uniform'; 
% Input.Marginals(counter).Parameters = [-1.5 1.5];
% Input.Marginals(counter).Bounds = [-1.5 1.5];



end



function Input  = NewMexico_calibrate()

% initialize counter of marginals
counter = 0;

% define marginals of yaw model
counter = counter+1;
AM11     = 0.445;
AM11_Std = 0.2;  % Standard deviation
% AM11_LB = 0; % Lower bound of truncated Gaussian distribution
% AM11_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM11';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM11, AM11_Std*abs(AM11)];
% Input.Marginals(counter).Bounds = [AM11_LB AM11_UB]; 

counter = counter+1;
AM12     = -1.78;
AM12_Std = 0.2;  % Standard deviation
% AM12_LB = 0; % Lower bound of truncated Gaussian distribution
% AM12_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM12';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM12, AM12_Std*abs(AM12)];
% Input.Marginals(counter).Bounds = [AM12_LB AM12_UB]; 

counter = counter+1;
AM13     = 1.63;
AM13_Std = 0.2;  % Standard deviation
% AM13_LB = 0; % Lower bound of truncated Gaussian distribution
% AM13_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM13';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM13, AM13_Std*abs(AM13)];

counter = counter+1;
AM14     = -0.0543;
AM14_Std = 0.2;  % Standard deviation
% AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM14';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM14, AM14_Std*abs(AM14)];

counter = counter+1;
AM15     = 0.367;
AM15_Std = 0.2;  % Standard deviation
% AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM15';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM15, AM15_Std*abs(AM15)];

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



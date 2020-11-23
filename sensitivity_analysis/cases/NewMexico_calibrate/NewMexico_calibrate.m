function Input  = NewMexico_calibrate()

% initialize counter of marginals
counter = 0;

% define marginals of Snel 3D Correction model
counter = counter+1;
SNEL3DFACTOR     = 3.1;
SNEL3DFACTOR_Std = 0.2;  % Standard deviation
SNEL3DFACTOR_LB = 1; % Lower bound of truncated Gaussian distribution
SNEL3DFACTOR_UB = 5;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'SNEL3DFACTOR';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [SNEL3DFACTOR, SNEL3DFACTOR_Std*abs(SNEL3DFACTOR)];
Input.Marginals(counter).Bounds = [SNEL3DFACTOR_LB SNEL3DFACTOR_UB]; 

counter = counter+1;
SNEL3DPOWER     = 2;
SNEL3DPOWER_Std = 0.2;  % Standard deviation
SNEL3DPOWER_LB = 1; % Lower bound of truncated Gaussian distribution
SNEL3DPOWER_UB = 3;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'SNEL3DPOWER';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [SNEL3DPOWER, SNEL3DPOWER_Std*abs(SNEL3DPOWER)];
Input.Marginals(counter).Bounds = [SNEL3DPOWER_LB SNEL3DPOWER_UB]; 


% % define marginals of yaw model
% counter = counter+1;
% AM11     = 0.445;
% AM11_Std = 0.2;  % Standard deviation
% % AM11_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM11_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM11';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM11, AM11_Std*abs(AM11)];
% % Input.Marginals(counter).Bounds = [AM11_LB AM11_UB]; 
% 
% counter = counter+1;
% AM12     = -1.78;
% AM12_Std = 0.2;  % Standard deviation
% % AM12_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM12_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM12';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM12, AM12_Std*abs(AM12)];
% % Input.Marginals(counter).Bounds = [AM12_LB AM12_UB]; 
% 
% % define marginals 
% % counter = counter+1;
% % Input.Marginals(counter).Name = 'CL';
% % Input.Marginals(counter).Airfoil = 'risoA121.dat';
% % Input.Marginals(counter).AlphaPert = [-10 30];
% % Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
% % Input.Marginals(counter).Type = 'Uniform'; 
% % Input.Marginals(counter).Parameters = [-1.5 1.5];
% % Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% % counter = counter+1;
% % Input.Marginals(counter).Name = 'CL';
% % Input.Marginals(counter).Airfoil = 'du91w250.dat';
% % Input.Marginals(counter).AlphaPert = [-10 30];
% % Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
% % Input.Marginals(counter).Type = 'Uniform'; 
% % Input.Marginals(counter).Parameters = [-1.5 1.5];
% % Input.Marginals(counter).Bounds = [-1.5 1.5];
% 
% % counter = counter+1;
% % Input.Marginals(counter).Name = 'CL';
% % Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
% % Input.Marginals(counter).AlphaPert = [-10 30];
% % Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
% % Input.Marginals(counter).Type = 'Uniform'; 
% % Input.Marginals(counter).Parameters = [-1.5 1.5];
% % Input.Marginals(counter).Bounds = [-1.5 1.5];% 
% 
% % counter = counter+1;
% % Input.Marginals(counter).Name = 'CD';
% % Input.Marginals(counter).Airfoil = 'naca64418_clean_Re0.7M.dat';
% % Input.Marginals(counter).AlphaPert = [0 20];
% % Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
% % Input.Marginals(counter).Type = 'Uniform'; 
% % Input.Marginals(counter).Parameters = [-1.5 1.5];
% % Input.Marginals(counter).Bounds = [-1.5 1.5];



end



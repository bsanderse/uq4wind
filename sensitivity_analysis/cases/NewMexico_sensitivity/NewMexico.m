function Input  = NewMexico()

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
Input.Marginals(counter).Parameters = [AM11, AM11_Std];
% Input.Marginals(counter).Bounds = [AM11_LB AM11_UB]; 

counter = counter+1;
AM12     = -1.78;
AM12_Std = 0.2;  % Standard deviation
% AM12_LB = 0; % Lower bound of truncated Gaussian distribution
% AM12_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM12';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM12, AM12_Std];
% Input.Marginals(counter).Bounds = [AM12_LB AM12_UB]; 

counter = counter+1;
AM13     = 1.63;
AM13_Std = 0.2;  % Standard deviation
% AM13_LB = 0; % Lower bound of truncated Gaussian distribution
% AM13_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM13';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM13, AM13_Std];
% Input.Marginals(counter).Bounds = [AM13_LB AM13_UB]; 

counter = counter+1;
AM14     = -0.0543;
AM14_Std = 0.2;  % Standard deviation
% AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM14';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM14, AM14_Std];
% Input.Marginals(counter).Bounds = [AM14_LB AM14_UB]; 

end



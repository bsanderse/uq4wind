function Input  = NewMexico()

% initialize counter of marginals
counter = 0;

% define marginals of yaw model
% counter = counter+1;
% AM11     = 0.445;
% AM11_Std = 1;  % Standard deviation
% % AM11_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM11_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM11';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM11, AM11_Std*abs(AM11)];
% % Input.Marginals(counter).Bounds = [AM11_LB AM11_UB]; 

% counter = counter+1;
% AM12     = -1.78;
% AM12_Std = 1*abs(AM12);  % Standard deviation
% % AM12_LB  = AM12 - 3*AM12_Std; % Lower bound of truncated Gaussian distribution
% % AM12_UB  = AM12 + 3*AM12_Std;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM12';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM12 AM12_Std];
% % Input.Marginals(counter).Bounds = [AM12_LB AM12_UB]; 
% 
% counter = counter+1;
% AM13     = 1.63;
% AM13_Std = 1;  % Standard deviation
% % AM13_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM13_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM13';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM13, AM13_Std*abs(AM13)];
% % Input.Marginals(counter).Bounds = [AM13_LB AM13_UB]; 
% 
% counter = counter+1;
% AM14     = -0.0543;
% AM14_Std = 1;  % Standard deviation
% % AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM14';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM14, AM14_Std*abs(AM14)];
% % Input.Marginals(counter).Bounds = [AM14_LB AM14_UB]; 

counter = counter+1;
AM15     = 0.367;
AM15_Std = 1;  % Standard deviation
% AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
Input.Marginals(counter).Name = 'AM15';
Input.Marginals(counter).Type = 'Gaussian'; 
Input.Marginals(counter).Parameters = [AM15, AM15_Std*abs(AM15)];
% Input.Marginals(counter).Bounds = [AM14_LB AM14_UB]; 


% counter = counter+1;
% AM21     = 0.0523;
% AM21_Std = 1;  % Standard deviation
% % AM11_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM11_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM21';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM21, AM21_Std*abs(AM21)];
% % Input.Marginals(counter).Bounds = [AM11_LB AM11_UB]; 
% 
% counter = counter+1;
% AM22     = -0.284;
% AM22_Std = 1;  % Standard deviation
% % AM12_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM12_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM22';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM22, AM22_Std*abs(AM22)];
% % Input.Marginals(counter).Bounds = [AM12_LB AM12_UB]; 
% 
% counter = counter+1;
% AM23     = 0.327;
% AM23_Std = 1;  % Standard deviation
% % AM13_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM13_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM23';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM23, AM23_Std*abs(AM23)];
% % Input.Marginals(counter).Bounds = [AM13_LB AM13_UB]; 
% 
% counter = counter+1;
% AM24     = -0.0134;
% AM24_Std = 1;  % Standard deviation
% % AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM24';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM24, AM24_Std*abs(AM24)];
% % Input.Marginals(counter).Bounds = [AM14_LB AM14_UB]; 
% 
% counter = counter+1;
% AM25     = 0.144;
% AM25_Std = 1;  % Standard deviation
% % AM14_LB = 0; % Lower bound of truncated Gaussian distribution
% % AM14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'AM25';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [AM25, AM25_Std*abs(AM25)];
% % Input.Marginals(counter).Bounds = [AM14_LB AM14_UB]; 
% 
% counter = counter+1;
% PH11     = -51.2;
% PH11_Std = 120;  % Standard deviation
% % PH11_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH11_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH11';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH11, PH11_Std];
% % Input.Marginals(counter).Bounds = [PH11_LB PH11_UB]; 
% 
% counter = counter+1;
% PH12     = 1009;
% PH12_Std = 120;  % Standard deviation
% % PH12_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH12_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH12';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH12, PH12_Std];
% %Input.Marginals(counter).Bounds = [PH12_LB PH12_UB]; 
% 
% counter = counter+1;
% PH13     = -1383;
% PH13_Std = 120;  % Standard deviation
% % PH13_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH13_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH13';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH13, PH13_Std];
% %Input.Marginals(counter).Bounds = [PH13_LB PH13_UB]; 
% 
% counter = counter+1;
% PH14     = 387;
% PH14_Std = 120;  % Standard deviation
% % PH14_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH14';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH14, PH14_Std];
% %Input.Marginals(counter).Bounds = [PH14_LB PH14_UB]; 
% 
% counter = counter+1;
% PH15     = -260;
% PH15_Std = 120;  % Standard deviation
% % PH14_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH15';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH15, PH15_Std];
% %Input.Marginals(counter).Bounds = [PH14_LB PH14_UB]; 
% 
% %  
% counter = counter+1;
% PH21     = 296;
% PH21_Std = 120;  % Standard deviation
% % PH11_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH11_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH21';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH21, PH21_Std];
% % Input.Marginals(counter).Bounds = [PH11_LB PH11_UB]; 
% 
% counter = counter+1;
% PH22     = 60.9;
% PH22_Std = 120;  % Standard deviation
% % PH12_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH12_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH22';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH22, PH22_Std];
% % Input.Marginals(counter).Bounds = [PH12_LB PH12_UB]; 
% 
% counter = counter+1;
% PH23     = -71.3;
% PH23_Std = 120;  % Standard deviation
% % PH13_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH13_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH23';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH23, PH23_Std];
% % Input.Marginals(counter).Bounds = [PH13_LB PH13_UB]; 
% 
% counter = counter+1;
% PH24     = -335;
% PH24_Std = 120;  % Standard deviation
% % PH14_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH24';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH24, PH24_Std];
% % Input.Marginals(counter).Bounds = [PH14_LB PH14_UB]; 
% 
% counter = counter+1;
% PH25     = 243;
% PH25_Std = 120;  % Standard deviation
% % PH14_LB = 0; % Lower bound of truncated Gaussian distribution
% % PH14_UB = 1;  % Upper bound of truncated Gaussian distribution 
% Input.Marginals(counter).Name = 'PH25';
% Input.Marginals(counter).Type = 'Gaussian'; 
% Input.Marginals(counter).Parameters = [PH25, PH25_Std];
% % Input.Marginals(counter).Bounds = [PH14_LB PH14_UB]; 

end
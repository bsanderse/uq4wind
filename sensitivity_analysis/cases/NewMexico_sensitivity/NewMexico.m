function Input  = NewMexico()

%% set type of distribution to describe the parameters that will be used to switch
    distr = 'Uniform';
% initialize counter of marginals
counter = 0;
%% define marginals of Snel 3D Correction model in case of performing 3D correction outside of Aeromodule

    switch distr

        case 'Gaussian'

        counter = counter+1;
        factor3D     = 2.2;
        factor3D_Std = 0.1;  % Standard deviation
        % factor3D_LB = 1.5; % Lower bound of truncated Gaussian distribution
        % factor3D_UB = 3.5;  % Upper bound of truncated Gaussian distribution 
        Input.Marginals(counter).Name = 'factor3D';
        Input.Marginals(counter).Type = distr; 
        Input.Marginals(counter).Parameters = [factor3D , factor3D_Std*abs(factor3D)];
        % Input.Marginals(counter).Bounds = [factor3D_LB factor3D_UB]; 

        counter = counter+1;
        exp3D     = 4;
        exp3D_Std = 0.1;  % Standard deviation
        % exp3D_LB = 2; % Lower bound of truncated Gaussian distribution
        % exp3D_UB = 5;  % Upper bound of truncated Gaussian distribution 
        Input.Marginals(counter).Name = 'exp3D';
        Input.Marginals(counter).Type = distr; 
        Input.Marginals(counter).Parameters = [exp3D, exp3D_Std*abs(exp3D)];
        % Input.Marginals(counter).Bounds = [exp3D_LB exp3D_UB]; 

        % %% define marginals of Snel 3D Correction model in case of including parameters in the specialist input
        % counter = counter+1;
        % SNEL3DFACTOR     = 3.1;
        % SNEL3DFACTOR_Std = 0.2;  % Standard deviation
        % SNEL3DFACTOR_LB = 1; % Lower bound of truncated Gaussian distribution
        % SNEL3DFACTOR_UB = 5;  % Upper bound of truncated Gaussian distribution 
        % Input.Marginals(counter).Name = 'SNEL3DFACTOR';
        % Input.Marginals(counter).Type = 'Gaussian'; 
        % Input.Marginals(counter).Parameters = [SNEL3DFACTOR, SNEL3DFACTOR_Std*abs(SNEL3DFACTOR)];
        % Input.Marginals(counter).Bounds = [SNEL3DFACTOR_LB SNEL3DFACTOR_UB]; 
        % 
        % counter = counter+1;
        % SNEL3DPOWER     = 2;
        % SNEL3DPOWER_Std = 0.2;  % Standard deviation
        % SNEL3DPOWER_LB = 1; % Lower bound of truncated Gaussian distribution
        % SNEL3DPOWER_UB = 3;  % Upper bound of truncated Gaussian distribution 
        % Input.Marginals(counter).Name = 'SNEL3DPOWER';
        % Input.Marginals(counter).Type = 'Gaussian'; 
        % Input.Marginals(counter).Parameters = [SNEL3DPOWER, SNEL3DPOWER_Std*abs(SNEL3DPOWER)];
        % Input.Marginals(counter).Bounds = [SNEL3DPOWER_LB SNEL3DPOWER_UB]; 

        case 'Uniform'

            counter = counter+1;

            factor3D_LB = 1.0; % Lower bound Uniform distribution
            factor3D_UB = 4;  % Upper bound of Uniform distribution 
            Input.Marginals(counter).Name = 'factor3D';
            Input.Marginals(counter).Type = distr; 
            Input.Marginals(counter).Bounds = [factor3D_LB factor3D_UB]; 

            counter = counter+1;

            exp3D_LB = 2; % Lower bound of Uniform distribution
            exp3D_UB = 5;  % Upper bound of Uniform distribution 
            Input.Marginals(counter).Name = 'exp3D';
            Input.Marginals(counter).Type = distr; 
            Input.Marginals(counter).Bounds = [exp3D_LB exp3D_UB]; 
    end
end



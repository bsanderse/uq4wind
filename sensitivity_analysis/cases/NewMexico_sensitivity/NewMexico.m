function Input  = NewMexico()

%% set type of distribution to describe the parameters that will be used to switch

%FixedParameters = P.FixedParameters;

% type = FixedParameters.correction;

% Select Correction type, type = 1 -->Snel, type = 2 --> Chaviaropoulos - Hansen
%or type = 3 for turbulent wake model sensitivity, Wilson model

type = 3 ; 

counter = 0;
    switch type
    
    % initialize counter of marginals
    

    %% If correction = 1 --> Snel
    case {1}
    
        counter = counter+1;
        factor3D     = 3;
        factor3D_Std = 0.1;  % Standard deviation
        % factor3D_LB = 1.5; % Lower bound of truncated Gaussian distribution
        % factor3D_UB = 3.5;  % Upper bound of truncated Gaussian distribution 
        Input.Marginals(counter).Name = 'factor3D';
        Input.Marginals(counter).Type = 'Gaussian'; 
        Input.Marginals(counter).Parameters = [factor3D , factor3D_Std*abs(factor3D)];
        % Input.Marginals(counter).Bounds = [factor3D_LB factor3D_UB]; 

        counter = counter+1;
        exp3D     = 2;
        exp3D_Std = 0.1;  % Standard deviation
        % exp3D_LB = 2; % Lower bound of truncated Gaussian distribution
        % exp3D_UB = 5;  % Upper bound of truncated Gaussian distribution 
        Input.Marginals(counter).Name = 'exp3D';
        Input.Marginals(counter).Type = 'Gaussian'; 
        Input.Marginals(counter).Parameters = [exp3D, exp3D_Std*abs(exp3D)];
        % Input.Marginals(counter).Bounds = [exp3D_LB exp3D_UB]; 

        %% If correction = 2 --> Chaviaropoulos - Hansen
    case {2}
        
        %% Chav. - Hansen model
        counter = counter+1;
        factor3D     = 2.2;
        factor3D_Std = 0.1;  % Standard deviation
        % factor3D_LB = 1.5; % Lower bound of truncated Gaussian distribution
        % factor3D_UB = 3.5;  % Upper bound of truncated Gaussian distribution 
        Input.Marginals(counter).Name = 'factor3D';
        Input.Marginals(counter).Type = 'Gaussian'; 
        Input.Marginals(counter).Parameters = [factor3D , factor3D_Std*abs(factor3D)];
        % Input.Marginals(counter).Bounds = [factor3D_LB factor3D_UB]; 

        counter = counter+1;
        exp3D     = 4;
        exp3D_Std = 0.1;  % Standard deviation
        % exp3D_LB = 2; % Lower bound of truncated Gaussian distribution
        % exp3D_UB = 5;  % Upper bound of truncated Gaussian distribution 
        Input.Marginals(counter).Name = 'exp3D';
        Input.Marginals(counter).Type = 'Gaussian'; 
        Input.Marginals(counter).Parameters = [exp3D, exp3D_Std*abs(exp3D)];
        % Input.Marginals(counter).Bounds = [exp3D_LB exp3D_UB];
        
    case {3}
        counter = counter +1;
        ATRANSITION = 0.38;
        ATRANSITION_std = 0.1;
        Input.Marginals(counter).Name = 'ATRANSITION';
        Input.Marginals(counter).Type = 'Gaussian'; 
        
        %mu = log(ATRANSITION^2/(sqrt(ATRANSITION^2+(ATRANSITION_std*abs(ATRANSITION))^2)));
        %std = sqrt(log(1+(ATRANSITION_std*abs(ATRANSITION))^2/ATRANSITION^2));
        
        Input.Marginals(counter).Parameters = [ATRANSITION, ATRANSITION_std*abs(ATRANSITION)];
        
    otherwise
        
        error('Invalid sensitivity case, check uncertain parameters')
        
    end
    
end



function Input  = NM80()


counter = 0;

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section03.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5]; 
Input.Marginals(counter).Bounds = [-.5 .5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section05.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 2; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5];
Input.Marginals(counter).Bounds = [-.5 .5];

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section08.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5];
Input.Marginals(counter).Bounds = [-.5 .5];% 

counter = counter+1;
Input.Marginals(counter).Name = 'CL';
Input.Marginals(counter).Airfoil = 'section10.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5];
Input.Marginals(counter).Bounds = [-.5 .5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section03.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 1; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5]; 
Input.Marginals(counter).Bounds = [-.5 .5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section05.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 2; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5];
Input.Marginals(counter).Bounds = [-.5 .5];

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section08.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 3; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5];
Input.Marginals(counter).Bounds = [-.5 .5];% 

counter = counter+1;
Input.Marginals(counter).Name = 'CD';
Input.Marginals(counter).Airfoil = 'section10.dat';
Input.Marginals(counter).AlphaPert = [-10 50];
Input.Marginals(counter).AirfoilIndex = 4; % unique ID corresponding to airfoil
Input.Marginals(counter).Type = 'Uniform'; 
Input.Marginals(counter).Parameters = [-.5 .5];
Input.Marginals(counter).Bounds = [-.5 .5];


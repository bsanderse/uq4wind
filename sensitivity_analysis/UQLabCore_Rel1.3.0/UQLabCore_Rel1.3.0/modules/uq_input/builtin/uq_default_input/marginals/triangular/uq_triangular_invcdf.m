function X = uq_triangular_invcdf(F, parameters)
% UQ_TRIANGULAR_INVCDF(F, parameters) calculates the inverse Cumulative 
% Density Function values of CDF values F of samples X  that follow  
% a triangular distribution with parameters specified in the vector
%  'parameters'

a = parameters(1);
b = parameters(2);
c = parameters(3);

X = zeros(size(F));
Fm = (c-a)/(b-a);

%% set the InvCDF to the appropriate value below Fm
idx = F <= Fm;
X(idx) = a + sqrt(F(idx) * (b-a)*(c-a));

%% set the InvCDF to the appropriate value above Fm
idx =  F > Fm ;
X(idx) = b - sqrt( (1 - F(idx))*(b-a)*(b-c));

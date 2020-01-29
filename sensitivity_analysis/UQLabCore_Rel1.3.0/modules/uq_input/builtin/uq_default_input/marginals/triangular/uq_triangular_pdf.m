function f = uq_triangular_pdf(X, parameters)
% UQ_TRIANGULAR_PDF(X, parameters) calculates the Probability Density Function
% values of samples X  that follow a triangular distribution with parameters 
% specified in the vector 'parameters'

a = parameters(1);
b = parameters(2);
c = parameters(3);

f = zeros(size(X));

%% set the PDF to 0 below the lower bound 
idx = X <= a;
f(idx) = 0;

%% set the PDF to the appropriate value in the valid range
idx = X > a & X < c ;
f(idx) = 2* (X(idx) - a) /(b-a)/(c-a) ;
idx = X == c;
X(idx) = 2/(b-a);
idx = X > c & X < b ;
f(idx) = 2* (b - X(idx)) /(b-a)/(b-c) ;

%% set the CDF to 1 above the bound 
idx = X >= b ;
f(idx) = 0;
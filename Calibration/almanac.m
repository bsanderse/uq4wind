function Y = almanac(Q,X)
% B: uncertain parameter, size p*1
% X: contains independent variable, size N*p: 

% y = b_1 + b_2*x 
Y = Q(1) + (X/12)*Q(2) + ((X/12).^2)*Q(3);

% y = b_1 + b_2*x + b_3*x^2
% Y = B(1) + X*B(2) + (X.^2)*B(3);
function Y = linearmodel(B,X)
% B: uncertain parameter, size p*1, or p*M
% X: contains independent variable, size N*p: 
% p: dimension of parameter vector
% N: number of locations
% M: number of parameter vectors to be evaluated

% output Y: size N*M 

% y = b_1 + b_2*x 
Y = B(1,:) + cos(X)*B(2,:);

% y = b_1 + b_2*x + b_3*x^2
% Y = B(1) + X*B(2) + (X.^2)*B(3);

function Y = linearmodel(X,B)
% X: independent variable, size N*p
% B: parameter, size p*1

% y = b_1 + b_2*x 
Y = B(1) + X*B(2);

% y = b_1 + b_2*x + b_3*x^2
% Y = B(1) + X(:,2)*B(2) + X(:,3).^2*B(3);
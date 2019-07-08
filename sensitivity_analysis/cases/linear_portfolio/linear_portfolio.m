function Y = linear_portfolio(X,a)
% linear combination of input values X with vector A:
% y = sum X_i * a_i
% X matrix of size Nxd 
% a vector of length d

% computing the response value
Y(:,1) = X*a; 

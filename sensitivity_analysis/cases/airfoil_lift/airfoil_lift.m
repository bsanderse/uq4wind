function Y = airfoil_lift(X)
% AIRFOIL_LIFT express lift in terms of parameters
%
% L = 0.5*rho*CL*V^2*c = 0.5*X(:,1)*X(:,2)^2
% with rho = 1 and c=1

% the response value:
Y(:,1) = 0.5*X(:,1).*(X(:,2)^2);

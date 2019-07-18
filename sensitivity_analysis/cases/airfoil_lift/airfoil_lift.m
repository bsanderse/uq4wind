function Y = airfoil_lift(X,P)
% AIRFOIL_LIFT express lift in terms of parameters
%
% L = 0.5*rho*CL*V^2*c = 0.5*X(:,1)*X(:,2)^2

rho = P(1);
c   = P(2);

% the response value:
Y(:,1) = 0.5*rho*X(:,1).*(X(:,2)^2)*c;

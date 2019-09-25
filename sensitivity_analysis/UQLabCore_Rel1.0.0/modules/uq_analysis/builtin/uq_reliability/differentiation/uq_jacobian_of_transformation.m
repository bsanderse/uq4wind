function J = uq_jacobian_of_transformation(X, U, Marginals, Copula, Options)
% J = UQ_JACOBIAN_OF_TRANSFORMATION (X,U,Marginals, Copula, Options):
%     computes the Jacobian matrix for a given point X and its transformed 
%     point U, with forward finite differences.
%     
% See also: UQ_FORM, UQ_GENERALISOPTRANSFORM

% Dimension
M = length(Marginals); 

% Initialize the transformation function:
[StandardMarginals(1:M).Type] = deal('Gaussian') ;
[StandardMarginals(1:M).Parameters] = deal([0,1]) ;
StandardCopula.Type = 'Gaussian';
StandardCopula.Parameters = eye(M);

% Compute the Jacobian with finite differences
h = Options.Gradient.h;
H = transpose(eye(M)*h);
MovedU = bsxfun(@plus, U, H);
MovedX = uq_GeneralIsopTransform(MovedU, StandardMarginals, StandardCopula, Marginals, Copula);

% Consider the transformed point as the output of a vector valued function:
MovedX = MovedX';

% Prepare a matrix that substracts the points in the original coordinates:
StaticX = repmat(X', 1, M);

% At the end we do:
% (h(i) is a vector with all zeros and h in the i-th position)
%
%         [ f_1(U + h(1)) - f_1(U), ... , f_1(U + h(M)) - f_1(U)]
% J = 1/h*[      ...                                   ...      ]
%         [ f_M(U + h(1)) - f_M(U), ... , f_M(U + h(M)) - f_M(U)]
%
J = (MovedX - StaticX)/h;
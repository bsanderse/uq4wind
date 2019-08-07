function c = getControlPoints(x,S,t,n)
%% This routine computes control points c for B-splines that corresponds to sampled data.
% Input arguments:
% 'x' is locations of known values of curve (e.g. locations where values of chord or twist are sampled)
% 'S' known function values at x
% 't' is the knot vector
% 'n' is the NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.

% Output arguments:
% 'c' is the values of control points at 'x' locations (at the sampled locations)

t = [t(1)*ones(1,n-1) t t(end)*ones(1,n-1)]; % padded knot vector obtained by padding n-1 elements at front and end. 
j = 0: numel(t)- n-1; % Index of B-spline from 0 =< j < numel(t)-n

B = zeros(numel(j),numel(j)); % Bspline matrix at 

for i = 1:numel(j)
    [B(:,i),x] = bspline_basis(j(i),n,t,x);
end
c = B\S'; % control points derived from sampled data 

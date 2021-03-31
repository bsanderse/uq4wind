function [B,t_padded] = getNURBSBasisMatrix(x,t,n)
%% This routine computes the Basis matrix for B-splines that corresponds to sampled data.
% Input arguments:
% 'x' is locations of known values of curve. Should be normalized between [0,1]
% 't' is the knot vector
% 'n' is the NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.

% Output arguments:
% 'B' is the basis matrix

t_padded = [t(1)*ones(1,n-1) t' t(end)*ones(1,n-1)]; % padded knot vector obtained by padding n-1 elements at front and end. 
j = 0: numel(t_padded)- n-1; % Index of B-spline from 0 =< j < numel(t)-n

B = zeros(numel(x),numel(j)); % Bspline matrix at x

for i = 1:numel(j)
    % evaluate basis function for entire vector x
    % bspline_basis returns a row vector which is stored as column vector
    % in B
    [B(:,i),x] = bspline_basis(j(i),n,t_padded,x);
end

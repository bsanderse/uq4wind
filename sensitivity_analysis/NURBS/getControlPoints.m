function c = getControlPoints(B,S)
%% This routine computes control points c for B-splines that corresponds to sampled data.
% Input arguments:
% 'B' Basis function matrix, each column corresponds to a certain basis
% function, each row to a certain position
% 'S' known function values at x

% Output arguments:
% 'c' is the values of control points at 'x' locations (at the sampled locations)

c = B\S(:); % control points derived from sampled data 

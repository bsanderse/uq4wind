%% This routine computes the B-splines y at locations x given a knot vector t and NURBS order n 

close all
clearvars
clc

n = 4; % 2 for linear B-splines, 3 for Quadratic, so on. Polynomial degree of B-spline is n-1.
t = [0 1/3 2/3 1]; % Original knot vector

% Knot vector is padded with (n-1) elements at front by repeating first
% element and at end by repeating the last element (n-1) times
t = [t(1)*ones(1,n-1) t t(end)*ones(1,n-1)];  
j = 0: numel(t)-n-1; % Index of knot span from 0 = < j < numel(t)-n. This is also equal to number of Bsplines.

hold all;
for i=1:numel(j)
    [y,x] = bspline_basis(j(i),n,t); % y is the Bspline value at locations x
    plot(x,y)
end
hold off;
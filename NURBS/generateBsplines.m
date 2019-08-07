%% This routine computes the B-splines y at locations x given a knot vector t and NURBS order n 

close all
clear all
clc

t = [0 1/3 2/3 1]; % knot vector, results in numel(t)+1 B-splines
n = 4; % 2 for linear B-splines, 3 for Quadratic, so on. Polynomial degree of B-spline is n-1. 
t = [t(1)*ones(1,n-1) t t(end)*ones(1,n-1)]; % padded knot vector obtained by padding(n-1) elements at front and end. 
j = 0: numel(t)-n-1; % Index of B-spline from 0 =< j < numel(t)-n
hold all;
for i=1:numel(j)
    [y,x] = bspline_basis(j(i),n,t); % y is the Bspline value at locations x
    plot(x,y)
end
hold off;

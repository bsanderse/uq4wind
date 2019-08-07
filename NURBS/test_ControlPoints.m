% This routine computes control points c for B-splines that corresponds to sampled data.

close all
clear all
clc

x = [1.5 15.42 32.25 52.75 63]/63; % locations of known values of curve (e.g. locations where values of chord or twist are sampled)
S = [3.28 4.62 3.70 2.55 1.37]; % known values at x
t = [0 1/3 2/3 1]; % knot vector, the number of resulting basis function is numel(t)+1
n = 3; % NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.

c = getControlPoints(x,S,t,n); % control points for NURB curve
t = [t(1)*ones(1,n-1) t t(end)*ones(1,n-1)]; % padded knot vector obtained by padding n-1 elements at front and end. 
j = 0: numel(t)- n-1; % Index of B-spline from 0 =< j < numel(t)-n

plot(x,c,'marker','o','linewidth',2) % plot control points
hold on
plot(x,S,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points
xu  = linspace(x(1), x(end), 100); % discrete query points
Sxu = zeros(1,numel(xu)); % 'Sxu' is the function values of NURB curve interpolated at xu

for i = 1:numel(j)
    [y,xu] = bspline_basis(j(i),n,t,xu);
    Sxu = Sxu + c(i)*y;
end

plot(xu,Sxu,'linewidth',2,'color','r')
legend('control points','sampled data', 'NURB curve')
hold off
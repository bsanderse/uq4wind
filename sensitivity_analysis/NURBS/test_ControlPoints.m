% This routine computes control points 'c' for B-splines that corresponds
% to sampled data 'S'.

close all
clearvars
clc

x = [1.5 15.42 32.25 52.75 63]/63; % locations of known values of curve (e.g. locations where values of chord or twist are sampled)
S = [3.28 4.62 3.70 2.55 1.37]; % Sampled data at x
n = 3; % NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.
t0 = [0 1/3 2/3 1]; % Original knot vector

% construct NURBS
[Bref, t] = getNURBSBasisMatrix(x,t0,n); % get basis matrix
c = getControlPoints(Bref,S); % control points for NURB curve

% evaluate NURBS
xu  = linspace(x(1), x(end), 100); % discrete query points
Bu  = getNURBSBasisMatrix(xu,t0,n);
Sxu = evalNURBS(Bu,c);

plot(x,c,'marker','o','linewidth',2) % plot control points
hold on
plot(x,S,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points
plot(xu,Sxu,'linewidth',2)
legend('control points','sampled data', 'NURBS curve')
hold off
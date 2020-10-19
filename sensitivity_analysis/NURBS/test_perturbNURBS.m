% This routine computes perturbed samples obtained from perturbing control
% points of baseline curve

close all
clearvars
clc

x = [1.5 15.42 32.25 52.75 63]/63; % locations of known values of curve (e.g. locations where values of chord or twist are sampled)
S = [3.28 4.62 3.70 2.55 1.37]; % known function values at x
t0 = [0 1/3 2/3 1]; % Original knot vectors 
n = 3; % NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.

% construct NURBS
[Bref, t] = getNURBSBasisMatrix(x,t0,n); % get basis matrix
c = getControlPoints(Bref,S); % control points for NURB curve

% evaluate NURBS
xu  = linspace(x(1), x(end), 100); % discrete query points
Bu  = getNURBSBasisMatrix(xu,t0,n);
Sxu = evalNURBS(Bu,c);

% add perturbation
samples = 10;
pc      = 0.1*ones(numel(c),1); % Amount of perturbation for each control points, for example, 0.1 corresponds to plus-minus 5% of perturbation.
randVec = rand(numel(c),samples) - 0.5; % Uniform distribution between [-0.5,0.5]

c_pert = c.*(1+pc.*randVec); 
S_pert = evalNURBS(Bu,c_pert);

plot(x,c,'marker','o','linewidth',2) % plot control points
hold on
plot(x,S,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points
plot(xu,Sxu,'linewidth',2)
plot(xu,S_pert,'color','k','linestyle','--')
legend('control points','sampled data', 'NURBS curve','random samples')
hold off
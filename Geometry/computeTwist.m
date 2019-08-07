% This routine computes the random samples of Twist vector using the purturbed
% control points of NURB curve

close all
clear all
clc

addpath('../NURBS') % location of NURBS routine directory
samples = 20; % Number of samples of perturbed twist

bladeLength = 38.8; % To normalize the blade length bw [0,1]

% Points along the blade length according to ECNAERO module input file 
interpolationLocations = [0 2 4 6 8 10 12 14 16 18 20 22 ...
                          24 26 28 30 32 34 36 37 38 38.4 38.8]/bladeLength; 
 
% Reference values of twist obtained from ECNAERO module input file                      
referenceTwist = [0 5.37 6.69 7.9 9.11 10.19 9.39 7.16 5.45 4.34 3.5 2.86 ... 
                 2.31 1.77 1.28 0.9 0.55 0.23 0.03 0.02 0.93 2.32 6.13]; 

% Plot to check the Twist variation along the blade span. This can be used to 
% select the knot locations            
plot(interpolationLocations,referenceTwist,'linewidth',2) 
hold on

t0 = [0 0.051 0.154 0.257 0.309 0.412 0.515 0.721 0.963 1]; % Knot vector b/w [0,1], the number of resulting basis function is numel(t)+1
n = 2; % NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.
sampleLocations = [0 2 6 10 12 16 20 28 37 38.8]/bladeLength; % Normalized between [0,1], numel(sampleLocations) =  numel(t0) + 1
sampledValues = [0 5.37 7.9 10.19 9.39 5.45 3.5 1.28 0.02 6.13]; % Vales of Chord at sampled location

c = getControlPoints(sampleLocations,sampledValues,t0,n); % control points for NURB curve
t = [t0(1)*ones(1,n-1) t0 t0(end)*ones(1,n-1)]; % padded knot vector obtained by padding n-1 elements at front and end. 
j = 0: numel(t)- n-1; % Index of B-spline from 0 =< j < numel(t)-n

plot(sampleLocations,c,'marker','o','linewidth',2) % plot control points
plot(sampleLocations,sampledValues,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points

sampleTwist = zeros(1,numel(interpolationLocations)); % 'sampleTwist' is the function values of NURB curve interpolated at sampleLocations

for i = 1:numel(j)
    [y,interpolationLocations] = bspline_basis(j(i),n,t,interpolationLocations);
    sampleTwist = sampleTwist + c(i)*y;
end

plot(interpolationLocations,sampleTwist,'linewidth',2,'color','g')
pc = 0.1*ones(numel(j),1);
sampleTwist = perturbNURBS(t0,n,interpolationLocations,c, pc,samples);
plot(interpolationLocations,sampleTwist','color','k','linestyle','--','HandleVisibility','off')
plot(interpolationLocations,sampleTwist(1,:),'color','k','linestyle','--')
legend('reference Twist','control points','sampled data', 'NURB curve','random samples')
hold off
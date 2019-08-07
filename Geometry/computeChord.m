% This routine computes the random samples of chord vector using the purturbed
% control points of NURB curve

close all
clear all
clc

addpath('../NURBS') % location of NURBS routine directory
samples = 20; % Number of samples of perturbed chord
bladeLength = 38.8; % To normalize the blade length bw [0,1]

% Points along the blade length according to ECNAERO module input file 
interpolationLocations = [0 2 4 6 8 10 12 14 16 18 20 22 ...
                          24 26 28 30 32 34 36 37 38 38.4 38.8]/bladeLength; 
                      
% Reference values of twist obtained from ECNAERO module input file                      
referenceChord = [2.42 2.48 2.65 2.81 2.98 3.14 3.17 2.99 2.79 2.58 2.38 ...
                  2.21 2.06 1.92 1.8 1.68 1.55 1.41 1.18 0.98 0.62 0.48 0.07]; 

% Plot to check the Chord variation along the blade span. This can be used to 
% select the knot locations            
plot(interpolationLocations,referenceChord,'linewidth',2) 
hold on

t0 = [0 0.15 0.3 0.51 0.7 0.87 0.95 1]; % Knot vector b/w [0,1], the number of resulting basis function is numel(t)+1
n = 3; % NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.
sampleLocations = [0 4 10 16 20 30 34 37 38.8]/bladeLength; % Normalized between [0,1], numel(sampleLocations) =  numel(t0) + 1
sampledValues = [2.42 2.65 3.14 2.79 2.38 1.68 1.41 0.98 0.07]; % Vales of Chord at sampled location

c = getControlPoints(sampleLocations,sampledValues,t0,n); % control points for NURB curve
t = [t0(1)*ones(1,n-1) t0 t0(end)*ones(1,n-1)]; % padded knot vector obtained by padding n-1 elements at front and end. 
j = 0: numel(t)- n-1; % Index of B-spline from 0 =< j < numel(t)-n

plot(sampleLocations,c,'marker','o','linewidth',2) % plot control points
plot(sampleLocations,sampledValues,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points

sampleChord = zeros(1,numel(interpolationLocations)); % 'sampleChord' is the function values of NURB curve interpolated at sampleLocations

for i = 1:numel(j)
    [y,interpolationLocations] = bspline_basis(j(i),n,t,interpolationLocations);
    sampleChord = sampleChord + c(i)*y;
end

plot(interpolationLocations,sampleChord,'linewidth',2,'color','r')

pc = 0.1*ones(numel(j),1);
sampleChord = perturbNURBS(t0,n,interpolationLocations,c, pc,samples);
plot(interpolationLocations,sampleChord','color','k','linestyle','--','HandleVisibility','off')
plot(interpolationLocations,sampleChord(1,:),'color','k','linestyle','--')
legend('reference Chord','control points','sampled data', 'NURB curve','random samples')
hold off

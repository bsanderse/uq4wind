function  samplesCurve = computeCurves(samples,index, randVec, pc, plotSamples,...
                                      interpolationLocations,referenceCurve,n,sampledIndices)
% This routine computes the random samples of Chord vector using the purturbed
% control points of NURB curve

% Input arguments
% 'samples' Number of samples of perturbed curve
% 'randVec' a samples-by-numOfControlPoints matrix of random numbers
% 'pc' is vector with each element bw [0,1] representing fraction of perturbation for each control point from their baseline value.
% 'plotSamples' 0 for no plot, 1 (default) to plot the generated samples, the baseline curve and the control points
% 'interpolationLocations' Points along the x-axis 
% 'referenceChord' Reference values of curve obtained from ECNAERO module input file
% 't0' Knot vector b/w [0,1], the number of resulting basis function is numel(t0)+1
% 'n' NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.
% 'sampledLocations' Locations where the value of chord is sample. NOTE: numel(sampleLocations) =  numel(t0) + 1
% 'sampledValues' Vales of chord at sampled location
% 'pc' is vector with each element bw [0,1] representing fraction of perturbation for control points from their baseline value. The numel(pc) = numel(sampledLocations)  

% Output argument
% 'samplesCurve' are samples of curves generated by perturbing the baseline control points. Each row corresponds to one sample.
min_X = min(interpolationLocations);
if(min_X < 0)
    interpolationLocations = interpolationLocations - min_X;
end
interpolationLocations = interpolationLocations./max(interpolationLocations); % Normalized between [0,1]
sampledLocations = interpolationLocations(sampledIndices);
sampledValues = referenceCurve(sampledIndices);

pc_mod = zeros(1,numel(sampledLocations));
pc_mod(index) = pc;

randVec_mod = zeros(1,numel(sampledLocations));
randVec_mod(index) = randVec;

if(n==2)
    t0 = sampledLocations;
elseif(n==3)
    center = ceil(length(sampledLocations)/2);
    t0 = [sampledLocations(1:center-1) sampledLocations(center+1:end)];
else
    disp("Use n = 2 (linear polynomial) or 3(quadratic polynomial)")
end

c = getControlPoints(sampledLocations,sampledValues,t0,n); % control points for NURB curve
t = [t0(1)*ones(1,n-1) t0 t0(end)*ones(1,n-1)]; % padded knot vector obtained by padding n-1 elements at front and end. 
j = 0: numel(t)- n-1; % Index of B-spline from 0 = < j < numel(t)-n

samplesCurve = zeros(1,numel(interpolationLocations)); % 'sampleChord' is the function values of NURB curve interpolated at sampleLocations

for i = 1:numel(j)
    [y,interpolationLocations] = bspline_basis(j(i),n,t,interpolationLocations);
    samplesCurve = samplesCurve + c(i)*y;
end

if plotSamples == 1
    % Plot to check the Chord variation along the blade span. This can be used to 
    % select the knot locations     
    figure
    plot(interpolationLocations,referenceCurve,'linewidth',2) 
    hold on
    % plot(sampledLocations,c,'marker','o','linewidth',2) % plot control points
    plot(sampledLocations,sampledValues,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points
    plot(interpolationLocations,samplesCurve,'linewidth',2,'color','b')
end
samplesCurve = perturbNURBS(t0,n,interpolationLocations,c, pc_mod,samples,randVec_mod);

if plotSamples == 1
    plot(interpolationLocations,samplesCurve','color','k','linestyle','--','HandleVisibility','off')
    plot(interpolationLocations,samplesCurve(1,:),'color','k','linestyle','--')
    legend('reference Curve','control points','sampled data', 'NURB curve','random samples')
    hold off
end
return
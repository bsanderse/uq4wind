function pertCurve = computeCurves(samples, index, randVec, pc, plotSamples,...
                                      interpolationLocations, referenceCurve, n, sampledIndices)
% This routine computes the random samples of a curve  using the perturbed
% control points of NURBS curve

% Input arguments
% 'samples' Number of samples of perturbed curve
% 'index' is a vector of integers prescribing the index of control points we want
% to perturb
% 'randVec' a samples-by-numOfControlPoints matrix of random numbers
% 'pc' is vector with each element between [0,1] representing fraction of perturbation for each control point from their baseline value.
% 'plotSamples' 0 for no plot, 1 (default) to plot the generated samples, the baseline curve and the control points
% 'interpolationLocations': Points along the x-axis where the curve is
% evaluated
% 'referenceCurve': Reference values of curve obtained from ECNAERO module
% input file, should have same length af interpolationLocations
% 'n' NURBS order: 2 for linear B-splines, 3 for Quadratic, so on. The polynomial degree of B-spline is n-1.
% 'sampledIndices': indices of interpolationLocations and referenceCurve
% which are used

% Output argument
% 'samplesCurve' are samples of curves generated by perturbing the baseline control points. Each row corresponds to one sample.



% %Example: -------------------
% samples = 1;
% index = [3:7];
% randVec = rand(samples,length(index));
% pc = 0.2*zeros(size(index));% plus minus 10 percent
% plotSamples =1;
% interpolationLocations = [0 2 4 6 8 10 12 14 16 18 20 22 ...
% 24 26 28 30 32 34 36 37 38 38.4 38.8];
% referenceCurve = [2.42 2.48 2.65 2.81 2.98 3.14 3.17 2.99 2.79 2.58 2.38 ...
% 2.21 2.06 1.92 1.8 1.68 1.55 1.41 1.18 0.98 0.62 0.48 0.07];
% n=3; % polynomial order 
% sampledIndices = [1 3 6 9 11 16 18 20 23]; % Location where you want to
% sample the reference curve
% samplesCurve = computeCurves(samples,index, randVec, pc, plotSamples,...
% interpolationLocations,referenceCurve,n,sampledIndices)
% %-------------------------------


%% set up NURBS
% minimum of angle of attack
min_X = min(interpolationLocations);
% shift to have the minimum larger than 0
if(min_X < 0)
    interpolationLocations = interpolationLocations - min_X;
end
% Normalized between [0,1]
interpolationLocations = interpolationLocations./max(interpolationLocations);

% 'sampledLocations' Locations where the value of chord is sampled. NOTE: numel(sampleLocations) =  numel(t0) + 1
% often this is the same as the interpolationlocations
sampledLocations = interpolationLocations(sampledIndices);
% 'sampledValues' Vales of curve at sampled location
sampledValues = referenceCurve(sampledIndices);

% construct knot vector
% 't0' Knot vector between [0,1], the number of resulting basis function is numel(t0)+1
if (n==2)
    t0 = sampledLocations;
elseif (n==3)
    center = ceil(length(sampledLocations)/2);
    t0 = [sampledLocations(1:center-1) sampledLocations(center+1:end)];
else
    error("Use n = 2 (linear polynomial) or 3 (quadratic polynomial)")
end

% get NURBS basis function matrix
[Bref, t] = getNURBSBasisMatrix(sampledLocations,t0,n); % get basis matrix
% get control points by solving the NURBS equation system
c = getControlPoints(Bref,sampledValues); % control points for NURBS curve

% now the NURBS curve is fully defined and can be evaluated at different
% positions
% 'samplesCurve' is the function values of NURBS curve interpolated at interpolationLocations
Bu  = getNURBSBasisMatrix(interpolationLocations,t0,n);
samplesCurve = evalNURBS(Bu,c);


%% now create perturbations
% create vector for perturbations
pc_mod = zeros(numel(sampledLocations),1);
% magnitude of perturbation that is used to scale the random numbers
pc_mod(index(:)) = pc;

% value of random variable
randVec_mod = zeros(numel(sampledLocations),1);
randVec_mod(index(:)) = randVec;

% this is a key step, where the perturbation is added to the baseline
c_pert    = c.*(1+pc_mod.*randVec_mod);
pertCurve = evalNURBS(Bu,c_pert);
% pertCurve = perturbNURBS(t0,n,interpolationLocations,c, pc_mod,samples,randVec_mod);


%% make plots
if plotSamples == 1
    %figure
    plot(interpolationLocations,referenceCurve,'linewidth',2) 
    hold on
    %plot(sampledLocations,c,'marker','o','linewidth',2) % plot control points
    plot(sampledLocations,sampledValues,'marker','x','markersize',8,'linestyle','none','linewidth',2) % plot sampled points
    plot(interpolationLocations,samplesCurve,'linewidth',2)
    plot(interpolationLocations,pertCurve,'linestyle','--')
    legend('reference curve','sampled data','NURBS approximation to reference','perturbed NURBS')
    hold off
end
return
function uq_LRA_display(LRAModel, outArray, varargin)
% uq_LRA_display(LRAModel)
% 
% This is implemented in the same manner as uq_PCE_display.
%
% Plot the log of the basis coefficients as a function of the number in the
% enumerated Psi_basis multiplied by the normalization factors for
% different ranks.
% 
% Colors and size of the circles correspond to the total degree of the
% corresponding polynomials

if ~LRAModel.Internal.Runtime.isCalculated
    fprintf('LRA module %s is not yet initialized!\nGiven Configuration Options:', LRAModel.Name);
    LRAModel.Options
    return;
end


%% Consistency checks and command line parsing
if ~exist('outArray', 'var')
    outArray = 1;
    if length(LRAModel.LRA) > 1
        warning('The selected LRA has more than 1 output. Only the 1st output will be printed');
        fprintf('You can specify the outputs you want to be displayed with the syntax:\n')
        fprintf('uq_display(LRAModule, OUTARRAY)\nwhere OUTARRAY is the index of desired outputs, e.g. 1:3 for the first three\n\n')
    end
end
if max(outArray) > length(LRAModel.LRA)
    error('Requested output range is too large') ;
end

%% Simply plot the true vs reconstructed experimental design
for ii = 1:length(outArray)
    oo = outArray(ii);
    X = LRAModel.ExpDesign.X;
    Y = LRAModel.ExpDesign.Y;
    YLRA = uq_evalModel(LRAModel, X);
    uq_figure('name', sprintf('Output #%i', oo))
    uq_plot(Y(:,oo),YLRA(:,oo),'+')
    axis equal tight
    hold on
    % plot a diagonal dashed line on top of the scatter plot
    XL = xlim;
    YL = ylim;
    minXY = min([XL(:); YL(:)]);
    maxXY = max([XL(:); YL(:)]);
    uq_plot([minXY maxXY],[minXY maxXY], '--','Color', [0.5 0.5 0.5])
    % Set title and label
    title(sprintf('LRA: ED fit for output %d',oo))
    xlabel('$\mathrm{Y_{ED}}$')
    ylabel('$\mathrm{Y_{LRA}}$')
end

function results = uq_reliability(CurrentAnalysis)
% UQ_RELIABILITY provides the starting point for any structural reliability
%     analysis in UQLab. It directs to the corresponding file for the actual
%     analysis
% 
% See also: UQ_MC, UQ_FORM, UQ_SORM, UQ_SR_IMPORTANCE_SAMPLING,
% UQ_SUBSETSIM, UQ_AKMCS

%% Set the options
Options = CurrentAnalysis.Internal;

%% Select the Method
switch lower(Options.Method)
    
    % Crude Monte Carlo Simulation
    case {'mc'} 
        
        if Options.Display > 0
            fprintf('\nStarting Crude Monte Carlo Analysis...\n');
        end
        
        results = uq_mc(CurrentAnalysis);
    
    % First Order Reliability Method    
    case 'form' 
        if Options.Display > 0
            fprintf('\nStarting FORM Analysis...\n');
        end
        
        results = uq_form(CurrentAnalysis);
    
    % Second Order Reliability Method
    case 'sorm' 
        if Options.Display > 0
            fprintf('\nStarting FORM/SORM Analysis...\n');
        end
        
        results = uq_sorm(CurrentAnalysis);
    
    % Importance Sampling    
    case {'is'}
        if Options.Display > 0
            fprintf('\nStarting Importance Sampling Analysis...\n');
        end
        
         results = uq_sr_importance_sampling(CurrentAnalysis);
        
    % Subset Simulation
    case 'subset'
        if Options.Display > 0
            fprintf('\nStarting Subset Simulation Analysis...\n');
        end
        
        results = uq_subsetsim(CurrentAnalysis);
        
    %AK-MCS
    case 'akmcs'
        if Options.Display > 0
            fprintf('\nStarting AK-MCS Analysis...\n');
        end
        
        results = uq_akmcs(CurrentAnalysis);
        
    otherwise
        error('The selected method "%s" does not exist', Options.Method);
        
end
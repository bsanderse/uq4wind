function [g, newMeta_flag] = uq_coupled_cmaesnonlconwrapper( d, current_analysis, update_flag, iteration )
% This flag will be updated to true if there is a new metamodel
newMeta_flag = false ;
switch lower( current_analysis.Internal.Method )
    case {'two level', 'two-level', 'ria', 'pma','qmc'}
        g = uq_twolevel_evalConstraints( d, current_analysis ) ;
    case {'decoupled', 'sora'}
        g = uq_sora_evalShiftedConstraint( d, current_analysis ) ;
    case 'sla'
        g = uq_sla_evalConstraints( d, current_analysis )' ;
    case 'deterministic'
        g = uq_deterministic_evalConstraints(d, current_analysis)' ;
end

if isfield(current_analysis.Internal.Constraints,'SoftConstModel')
    gs = uq_evalSoftConstraint(d, current_analysis )' ;
else
    % Give a negative dummy value just to pass the test on update flag
    gs = -1 ;
end

% update_flag is true if the current point improves the current best
% point
% Enrichment is made is update_flag is true and the current design is
% feasible
% To improve either add a relaxed criterion on g (e.g. g<0.1*Pf_threshold or some sigmoid function)
% Use restart each time the enrichment is succesfull
try_enr_unfeasible = false ;
switch current_analysis.Internal.Constraints.Enrichment.LocalTryEnr
    case 1
        try_enr = true ;
    case 2
        try_enr = (update_flag == 1) ;
    case 3
        try_enr = ( update_flag == 1 && all(gs<=0) && all(g<=0) );
    case 4
        %                 try_enr = ( update_flag == 1 && all(gs<=0) && all(g<=0) ) ;
        try_enr = ( update_flag == 1 && all(gs<=0) ) ;
        try_enr_unfeasible = true ;
end

if try_enr
    for ii = 1:current_analysis.Internal.Constraints.Enrichment.LocalMaxAdded
        current_analysis.Internal.Runtime.Model = current_analysis.Internal.Constraints.Model ;
        TH = 0.1 ;
        % For all applications
%         if iteration <= 150
%             TH = 1 ;
%         elseif iteration >= 150 && iteration < 300
%             TH = 0.5 ;
%         elseif iteration >= 300 && iteration < 475
%             TH = 0.1 ;
%         elseif iteration >= 475 && iteration < 500
%             TH = 0.05 ;
%         else
%             TH = Inf ;
%         end
        
        % For bracket structure
%         if iteration <= 150
%             TH = 10 ;
%         elseif iteration >= 150 && iteration < 300
%             TH = 5 ;
%         elseif iteration >= 300 && iteration < 475
%             TH = 2 ;
%         elseif iteration >= 475 && iteration < 500
%             TH = 1 ;
%         else
%             TH = Inf ;
%         end
        
                
        % For dome application - as
%         if iteration <= 300
%             TH = 0.5 ;
%         elseif iteration >= 300 && iteration < 500
%             TH = 0.25 ;
%         elseif iteration >= 500 && iteration < 700
%             TH = 0.1 ;
%         elseif iteration >= 700 && iteration < 1000
%             TH = 0.05 ;
%         else
%             TH = Inf ;
%         end
%         std_gX = current_analysis.Internal.Runtime.std_gX ;
        if ~isfield(current_analysis.Internal.Runtime.ActiveMeta,'sigma_C')
            xcandidate = current_analysis.Internal.Runtime.ActiveMeta.xcandidate ;
            ycandidate = uq_evalModel(current_analysis.Internal.Constraints.Model, xcandidate);
            std_C = std(ycandidate) ;
            current_analysis.Internal.Runtime.ActiveMeta.sigma_C = std_C ;
        end
        % Minus case: g - 1.96 sigma
        minus_ca = current_analysis ;
        MMinus.mHandle = @(X) uq_g_minus(current_analysis.Internal.Runtime.Model,X);
        MMinus.isVectorized = true;
        minus_ca.Internal.Constraints.Model = uq_createModel(MMinus, '-private') ;
        g_minus = uq_twolevel_evalConstraints(d, minus_ca) ;
%          std_minus = current_analysis.Internal.Runtime.std_gX ;
       
        % Plus case: g + 1.96 sigma
        plus_ca = current_analysis ;
        MPlus.mHandle = @(X) uq_g_plus(current_analysis.Internal.Runtime.Model,X);
        MPlus.isVectorized = true;
        plus_ca.Internal.Constraints.Model = uq_createModel(MPlus, '-private') ;
        g_plus = uq_twolevel_evalConstraints(d, plus_ca) ;
%          std_plus = current_analysis.Internal.Runtime.std_gX ;
        
        % Get the actual constraints
        current_analysis.Internal.Runtime.isnewpoint = true ;
        Target = uq_recordconstraints(g, current_analysis) ;
        Target_minus = uq_recordconstraints(g_minus, current_analysis) ;
        Target_plus = uq_recordconstraints(g_plus, current_analysis) ;
        
%         Target = Target ./  std_gX;
%         Target_plus = Target_plus ./ std_plus ;
%         Target_minus = Target_minus ./ std_minus ;

        % Get the metamodel ?
        current_analysis.Internal.Constraints.Model = current_analysis.Internal.Runtime.Model ;
        
        switch lower(current_analysis.Internal.Reliability.Method)
            case {'mcs', 'subset','is'}
                switch lower(current_analysis.Internal.Optim.ConstraintType)
                    case 'pf'
                        crit =  abs(Target_minus - Target_plus)./(abs(Target) ) ;
                    case 'beta'
                        crit =  abs(Target_minus - Target_plus)./(abs(Target) + 1 ) ;
                end
            case 'qmc'
%                 crit =  abs(Target_minus - Target_plus)./(abs(Target) + 1) ;
crit = abs(Target_plus - Target_minus)./ current_analysis.Internal.Runtime.ActiveMeta.sigma_C ;
        end
        crit = mean(crit) ;
        % Record the constraint
        current_analysis.Internal.Runtime.LocalEnrichConv = ...
            [current_analysis.Internal.Runtime.LocalEnrichConv, [crit; iteration; TH]] ;
        if all(crit <= TH)
            break ;
        end
        
        if try_enr_unfeasible
            if all(g <= 0) && all(g_minus <= 0) && all(g_plus <= 0)
                % This means we are far enough from the limit-state surface so
                % no need to enrich even if the Beta bounds is larger than
                % threshold
                break;
            end
        end
        
        XC = current_analysis.Internal.Runtime.XC ;
        if iscell(XC)
            XC = vertcat(XC{:});
        end
        
        [M_X, M_var] = uq_evalModel(current_analysis.Internal.Constraints.Model, XC);
        M_s = sqrt(M_var) ;
        %Limit-state options
        LSOptions = current_analysis.Internal.LimitState ;
        TH = LSOptions.Threshold ;
        switch LSOptions.CompOp
            case {'<', '<=', 'leq'}
                g_X  = M_X  - repmat(TH,size(M_X,1),1);
            case {'>', '>=', 'geq'}
                g_X  = repmat(TH,size(M_X,1),1) - M_X;
        end
        lf = uq_LF_U(g_X,M_s);
        if size(lf,2) > 1
            [lf,indcons] = max(lf,[],2) ;
        end
        [~, lfidx] = max(lf);
        
        Xadd = XC(lfidx,:);
        Yadd = uq_evalModel(current_analysis.Internal.Constraints.FullModel, Xadd);
        if isnan(Yadd)
            % Take the second best point
            sorted = sortrows(Xadded, lf, size(Xadd,2)+1);
            Xadd = sorted(end-1,1:size(Xadd,2)) ;
            Yadd = uq_evalModel(current_analysis.Internal.Constraints.FullModel, Xadd);
        end
        % Now if it is still NaN, skip the enrichment
        if isnan(Yadd)
            Xadd = [] ; Yadd = [] ;
            g = - Inf ;
        else
            metaopts = current_analysis.Internal.Constraints.Model.Options ;
            
            metaopts.ExpDesign.X = [ metaopts.ExpDesign.X; Xadd] ;
            metaopts.ExpDesign.Y = [ metaopts.ExpDesign.Y; Yadd] ;
            current_analysis.Internal.Constraints.Model = uq_createModel(metaopts) ;
            newMeta_flag = true ;
            
            Nadded = size(metaopts.ExpDesign.X,1) - current_analysis.Internal.Runtime.ActiveMeta.Ntotal;
            fprintf(['Active Metamodel RBDO - Coupled: ',num2str(Nadded), ' samples added\n']) ;
            
            % Update the std in the augmented space...
            xcandidate = current_analysis.Internal.Runtime.ActiveMeta.xcandidate ;
            ycandidate = uq_evalModel(current_analysis.Internal.Constraints.Model, xcandidate);
            std_C = std(ycandidate) ;
            current_analysis.Internal.Runtime.ActiveMeta.sigma_C = std_C ;
            g = uq_twolevel_evalConstraints( d, current_analysis ) ;
            if any(g >0)
                % This means no need to enrich further (in the case we didn't reach yet the maximumn enrichment size) because the point is likely not
                % feasible.
                break;
            end
        end
    end
end

if isfield(current_analysis.Internal.Constraints,'SoftConstModel')
    g = [ g; gs] ;
end

% Save current model response in terms of Pf, beta and g (if deterministic)
current_analysis.Internal.Runtime.isnewpoint = true ;
RecordedConstraints = uq_recordconstraints( g, current_analysis ) ;
current_analysis.Internal.Runtime.RecordedConstraints = [current_analysis.Internal.Runtime.RecordedConstraints; RecordedConstraints ] ;

end

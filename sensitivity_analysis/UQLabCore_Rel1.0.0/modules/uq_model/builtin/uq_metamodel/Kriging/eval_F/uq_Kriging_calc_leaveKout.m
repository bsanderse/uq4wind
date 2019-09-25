function [ errors, additional_metrics, auxMatrices] = ...
    uq_Kriging_calc_leaveKout(Nclasses, N, Y, R, F, run_case )
%UQ_KRIGING_CALC_LEAVEKOUT: Returns the Leave-K-Out error along with
% some additional metrics if requested.
%
% See also UQ_KRIGING_CALCULATE, UQ_KRIGING_EVAL_J_OF_THETA_CV

errors = cell(Nclasses,1) ;
sigma_i_Sq = cell(Nclasses,1) ;

% Produce the indices of each part that is then going to be treated either
% as a training or a validation set
nElemsPerPart = floor(N/Nclasses) ;
nIndices = repmat(nElemsPerPart, 1 , Nclasses);
nIndices(end) = N - nElemsPerPart * (Nclasses-1) ;

% randomly permute Y (only when it's not K = 1 ) 
if Nclasses < N
    randInd = randperm(N);
else
    randInd = 1:N;
end

randInd = mat2cell(randInd,1, nIndices );

auxMatrices = uq_Kriging_calc_auxMatrices( R, F, Y, run_case );
        
FTRinv = auxMatrices.FTRinv ;
FTRinvF = auxMatrices.FTRinvF ;

% Try to calculate B1 matrix as efficiently as possible
L = auxMatrices.cholR;
if ~isnan(L)
    FTRinvF_inv = auxMatrices.FTRinvF_inv ;
    if ~isempty(FTRinvF_inv)
        MM = eye(N) - F * FTRinvF_inv * FTRinv  ;
        B1 = L \ (transpose(L) \  MM );
    else
        MM = eye(N) - F * (FTRinvF \ FTRinv)  ;
        B1 = L \ (L.' \  MM );
    end
else
    Rinv = auxMatrices.Rinv ;
    FTRinvF_inv = auxMatrices.FTRinvF_inv;
    if ~isempty(FTRinvF_inv)
        MM = eye(N) - F * FTRinvF_inv * FTRinv  ;
        B1 = Rinv * MM ;
    else
        MM = eye(N) - F * (FTRinvF \ FTRinv)  ;
        B1 = Rinv * MM ;
    end
end

LOOmean = zeros(Nclasses, 1);
LOOsd = zeros(Nclasses, 1);
fullInd = 1:N;


switch nargout
    case {0,1}
        % this is normally the case during optimization of the
        % Cross-Validation objective function. In that case only the
        % Leave-K-Out error is calculated and no other error metrics
        if Nclasses ~= N
            
            for jj = 1 : Nclasses
                % get the indices of the current training and validation set
                %indTrain = [randInd{1:jj-1} randInd{jj+1:K} ];
                indValidate = randInd{jj} ;
                
                tidx = true(size(fullInd)) ;
                tidx(indValidate) = false ;
                
                indTrain = fullInd(tidx) ;
                
                
                y = Y(indTrain ) ;
                y0 = Y(indValidate) ;
                
                % Calculate leave one (or K) out predictions in a computationally efficient
                % manner based on [Dubrule 1983]
                if numel(indValidate) == 1
                    muY0 = -(B1(indValidate, indValidate) \ B1(indValidate, indTrain)) * y ;
                else
                    muY0 = -(B1(indValidate, indValidate) \ B1(indValidate, indTrain)) * y ;
                end
                errors{jj} = (muY0 - y0).^2 ;
            end
        else
            % case of leave-one-out [Dubrule 1983]
            sigma2Y0 = 1./diag(B1);
            LOOmean = -(bsxfun(@times,sigma2Y0,(B1*Y)) -Y);
            errors = num2cell((LOOmean - Y).^2);
        end
        return
    case {2,3}
        % this is the case when various Leave-K-Out error metrics are
        % requested. Normally this is executed on an already trained
        % Kriging metamodel regardless of the the estimation method that
        % was selected
        
        % case of general leave-K-out
        if Nclasses ~= N
            for jj = 1 : Nclasses
                % get the indices of the current training and validation set
                indValidate = randInd{jj} ;
                
                tidx = true(size(fullInd)) ;
                tidx(indValidate) = false ;
                
                indTrain = fullInd(tidx) ;
                
                
                y = Y(indTrain ) ;
                y0 = Y(indValidate) ;
                
                % Calculate leave one (or K) out predictions in a computationally efficient
                % manner based on [Dubrule 1983]
                if numel(indValidate) == 1
                    sigma2Y0 = 1/ B1(indValidate, indValidate) ;
                    muY0 = -(B1(indValidate, indValidate) \ B1(indValidate, indTrain)) * y ;
                else
                    muY0 = -(B1(indValidate, indValidate) \ B1(indValidate, indTrain)) * y ;
                    sigma2Y0 = B1(indValidate, indValidate) \ eye( numel(indValidate))  ;
                end
                errors{jj} = (muY0 - y0).^2 ;
                sigma_i_Sq{jj} = diag(sigma2Y0) ;
                
                if numel(indValidate) == 1
                    LOOmean(jj) = muY0;
                    LOOsd(jj) = sqrt(sigma_i_Sq{jj});
                end
            end
            
        else % case of leave-one-out [Dubrule 1983]
            sigma2Y0 = 1./diag(B1);
            muY0 = -(bsxfun(@times,sigma2Y0,(B1*Y)) -Y);
            LOOmean = muY0;
            LOOsd = sqrt(sigma2Y0);
            sigma_i_Sq = num2cell(sigma2Y0);
            errors = num2cell((LOOmean - Y).^2);
        end
        additional_metrics.sigma_i_Sq = sigma_i_Sq;
        additional_metrics.LOOmean = LOOmean;
        additional_metrics.LOOsd = LOOsd;
        
    otherwise
        error('Too many output arguments were requested!')
end



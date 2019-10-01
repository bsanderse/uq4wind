function [alphas,exitflag,lambda,K] = uq_SVC_optimizer( X, Y, current_model )
%UQ_SVC_OPTIMIZER Solve the QPC problem
%   Detailed explanation goes here

% Obtain the current output
current_output = current_model.Internal.Runtime.current_output ;
% N = size(current_model.ExpDesign.Y, 1);
C = current_model.SVC(current_output).C ;
%<TMPTRANSPOSE>
X = X.';
Y = Y.';
N = length(Y);
% Set up the quadratic optimization problem
switch lower(current_model.Internal.SVC(current_output).Penalization)
    case 'linear'
        % Equality constraints
        Aeq = Y ;
        beq = 0 ;
        % Bounds
        lb = zeros(N,1);
        ub = C * ones(N,1) ;
%         vectorC = zeros(N,1);
        % Calculate the matrices of the quadratic objective function
        theta = current_model.SVC(current_output).Kernel.Params;
        K_Family = current_model.SVC(current_output).Kernel.Family;
        K = uq_SVC_eval_K( X, X, theta, K_Family);
        H = K .* (Y'*Y);
        if norm(H-H',inf) > eps    % Parfaitement symétriser H si nécessaire -- Asymétrie due aux approx. num.
            norm(H-H',inf);
            H = (H+H')/2;
        end
        % Add an option for ridge here - Never used such thing
        
        % Add option for unbalanced data
        f = -ones(N,1);
        
        % Add option for starting point
        %     x0 = given options or x0 = zeros(N,1);
        
    case 'quadratic'
        % Equality constraints
        Aeq = Y ;
        beq = 0 ;
        % Bounds
        lb = zeros(N,1);
        ub = Inf * ones(N,1) ;
        vectorC = 1/C * ones(N,1);
        % Calculate the matrices of the quadratic objective function
        theta = current_model.SVC(current_output).Kernel.Params;
        K_Family = current_model.SVC(current_output).Kernel.Family;
        K = uq_SVC_eval_K( X, X, theta, K_Family);
        H = K .* (Y'*Y);
        if norm(H-H',inf) > eps    % Parfaitement symétriser H si nécessaire -- Asymétrie due aux approx. num.
            norm(H-H',inf);
            H = (H+H')/2;
        end
        H = H + diag(vectorC);
        
        % Add an option for ridge here - Never used such thing
        
        % Add option for unbalanced data
        f = -ones(N,1);
        
        % Add option for starting point
        %     x0 = given options or x0 = zeros(N,1);
end
current_model.Internal.H = H;
current_model.Internal.f = f;
current_model.Internal.Aeq = Aeq;
current_model.Internal.beq = beq;

% Solve the optimization problem
switch lower (current_model.Internal.SVC(1).QPSolver)
    case 'ip'
        %     options = optimset('Algorithm','interior-point-convex','TolFun', 1e-10,'Display','none' );
        options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','none','TolX',1e-15,'TolCon',1e-15);
    case 'as'
        options = optimoptions('quadprog','Algorithm','active-set','Display','none','TolX',1e-15,'TolCon',1e-15);
    case 'smo'
        error('SMO Optimization method is not yet implemented \n')
    otherwise
        error('Unknown Optimization method \n')
end
% Run the optimization
[alphas,~,exitflag,~,lambda]= quadprog(H,f,[],[],Aeq,beq,lb,ub,[],options);
% Usiing QPC solver (Quadratic programming in C)
% [alphas,exitflag,lambda] = qpip(H,f,[],[],Aeq,beq,lb,ub,0,[],0);
% lambda.eqlin = lambda.equality;

end
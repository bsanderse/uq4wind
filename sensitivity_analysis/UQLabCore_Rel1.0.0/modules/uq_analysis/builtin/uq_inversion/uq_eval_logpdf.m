function value = uq_eval_logpdf(X,Input)
% UQ_EVAL_LOGPDF calculates the logarithm of the multivariate probability 
%     density function value that corresponds to samples of a random vector 
%     collected in X with marginal distributions and copula specified by 
%     the Input structure (or uq_default_input object) Input. Currently 
%     only Gaussian copulas are implemented.
%
% See also UQ_EVAL_PDF, UQ_ALL_PDF, UQ_ALL_CDF, UQ_ALL_INVCDF

marginals = Input.Marginals;
copula    = Input.Copula;
logXpdf   = log(uq_all_pdf(X,marginals));

switch lower(copula.Type)

    case 'independent'
        % if the components of X are independent then the value of the log
        % of the multivariate PDF is simply the sum of the univariate log PDFs
        value = sum(logXpdf,2);

    case 'gaussian'
        % the components of X have dependency that is described by a
        % Gaussian copula
        R = copula.Parameters;
        M = length(marginals);

        % Specify the U marginals in the standard normal space
        [U_marginals(1:M).Type]       = deal('Gaussian') ;
        [U_marginals(1:M).Parameters] = deal([0,1]) ;

        % Marginal CDFs
        Xcdf = uq_all_cdf(X,marginals);

        % Gaussian inverse CDF
        invXcdf = uq_all_invcdf(Xcdf,U_marginals);
        
        % Finally calculate the multivariate PDF value based on the 
        % logairthm of Sklar's theorem
        % log(f_joint(x)) = log(c(F(x))) + sum(log(f(x_i)))
        logDetR =  2 * sum(log(diag(chol(R))));
        Cgauss =  (1/2)*logDetR -1/2*dot((invXcdf*(inv(R) - eye(M))),invXcdf, 2);
        value = Cgauss + sum(logXpdf,2);

    otherwise 
        error('Copulas of type "%s" are not supported',copula.Type);
end
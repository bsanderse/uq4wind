function value = uq_eval_pdf(X,Input)
% value = UQ_EVAL_PDF(X,Input):
%     Calculates the multivariate probability density function value that 
%     corresponds to samples of a random vector collected in X with 
%     marginal distributions and copula specified by the Input structure 
%     (or uq_default_input object) Input. Currently only Gaussian copulas
%     are implemented.
%
% See also UQ_ALL_PDF, UQ_ALL_CDF, UQ_ALL_INVCDF

marginals = Input.Marginals;
copula    = Input.Copula;
Xpdf      = uq_all_pdf(X,marginals);

switch lower(copula.Type)

    case 'independent'
        % if the components of X are independent then the value of the
        % multivariate PDF is simply the product of the univariate PDFs
        value = prod(Xpdf,2);

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
        
        % Finally calculate the multivariate PDF value based on Sklar's
        % theorem
        % f_joint(x) = c(F(x))*prod(f(x_i))
        Cgauss =  1/det(R)^(1/2)*exp(-1/2*dot((invXcdf*(inv(R) - eye(M))),invXcdf, 2));
        value = Cgauss.*prod(Xpdf,2);

    otherwise 
        error('Copulas of type "%s" are not supported',copula.Type);
end
% UQ_CREATEMODEL  create a UQLab MODEL object. 
%    myModel = UQ_CREATEMODEL(MODELOPTS) creates and stores in the current
%    UQLab session an object of type MODEL based on the configuration
%    options given in the MODELOPTS structure. 
%
%    myModel = UQ_CREATEMODEL(MODELOPTS, '-private') creates an object of
%    type MODEL based on the configuration options given in the MODELOPTS
%    structure, without saving it in the UQLab session. This is useful,
%    e.g., in case of adaptive algorithms which create a large number of
%    temporary MODEL objects that do not need to be stored in memory.
%
%    For more details about the available MODEL types, please refer to
%    model-specific initialization: 
%       help <a href="matlab:uq_ModelOptions">uq_ModelOptions</a>      - Computational model (e.g. MATLAB functions)
%       help <a href="matlab:uq_PCEOptions">uq_PCEOptions</a>        - Polynomial Chaos Expansion metamodel options
%       help <a href="matlab:uq_KrigingOptions">uq_KrigingOptions</a>    - Kriging metamodel options
%       help <a href="matlab:uq_PCKOptions">uq_PCKOptions</a>        - Polynomial Chaos Kriging metamodel options
%       help <a href="matlab:uq_LRAOptions">uq_LRAOptions</a>        - Low Rank Approximations metamodel options
%       help <a href="matlab:uq_SVCOptions">uq_SVCOptions</a>        - Support Vector Classification options
%       help <a href="matlab:uq_SVROptions">uq_SVROptions</a>        - Support Vector Regression options
%       help <a href="matlab:uq_UQLinkOptions">uq_UQLinkOptions</a>     - UQLink options
%
%    See also: uq_ModelOptions, uq_evalModel, uq_getModel, uq_listModels,
%              uq_selectModel, uq_createInput, uq_createAnalysis, uq_doc
%
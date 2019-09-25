% UQ_CREATEANALYSIS   create a UQLab ANALYSIS object.
%    myAnalysis = UQ_CREATEANALYSIS(ANALYSISOPTS) creates and stores in the
%    current UQLab session an object of type ANALYSIS based on the configuration
%    options given in the ANALYSISOPTS structure. 
%
%    myAnalysis = UQ_CREATEANALYSIS(ANALYSISOPTS, '-private') creates an
%    object of type ANALYSIS based on the configuration options given in
%    the ANALYSISOPTS structure, without saving it in the UQLab session.
%    This is useful, e.g., in parametric studies which create a large
%    number of temporary ANALYSIS objects that do not need to be stored in
%    memory. 
%
%    For more details about the available model types, please refer to
%    analysis-specific initialization: 
%       help <a href="matlab:uq_ReliabilityOptions">uq_ReliabilityOptions</a> - Reliability analysis (rare events estimation)
%       help <a href="matlab:uq_SensitivityOptions">uq_SensitivityOptions</a> - Sensitivity analysis
%       help <a href="matlab:uq_InversionOptions">uq_InversionOptions</a>   - Bayesian calibration options

%
%    See also: uq_getAnalysis, uq_listAnalyses, uq_selectAnalysis,
%              uq_createInput, uq_createModel,  uq_doc
%
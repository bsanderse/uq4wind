function pass = uq_uqlink_test_Archiving( level )
% UQ_UQLINK_TEST_ARCHIVING tests that the archiving options of the module
% work fine
%
% See also: UQ_READ_SSBEAMDEFLECTION.m UQ_UQLINK_TEST_POSSIBLECASES.m

% parameters
pass = true ;
eps_th = 1e-5 ;
% X = [b h L E p]; Y = V = beam deflection
X = [0.15 0.3 5 30000e6 1e4] ;
 
Ytrue = 5 * X(:,5) .* X(:,3).^4 ./ (32 * X(:,4) .* X(:,1) .* X(:,2).^3) ;

if nargin < 1
    level = 'normal';
end
fprintf(['\nRunning: |' level '| uq_uqlink_test_Archiving...\n']);

uqlab('-nosplash')

%% CASE 1: SAVE IS DEFAULT. CHECK THAT THE ZIP EXIST
EXECBASENAME = 'uq_SimplySupportedBeam_v2';
if ispc
    EXECSUFFIX = 'win';
elseif isunix
    if ~ismac
        EXECSUFFIX = 'linux';
    else
        EXECSUFFIX = 'mac';
    end
end
EXECNAME = [EXECBASENAME '_' EXECSUFFIX];

Mopts.Type = 'UQLink' ;
Mopts.Name = 'my Beam' ;
Mopts.Command = [fullfile(uq_rootPath,'modules','uq_model','builtin','uq_uqlink','test','Case2', EXECNAME), ' SSB_Input_v2_1.inp', ' SSB_Input_v2_2.inp'] ;    
Mopts.ExecutionPath = fullfile(uq_rootPath,'modules','uq_model','builtin','uq_uqlink','test','Case2') ;
Mopts.Template = {'SSB_Input_v2_1.inp.tpl', 'SSB_Input_v2_2.inp.tpl'} ;
Mopts.Output.Parser= 'uq_read_SSBeamDeflection' ;
Mopts.Output.FileName = 'output.out' ;
Mopts.Display = 'quiet' ;

myModel = uq_createModel(Mopts) ;
Yval = uq_evalModel(myModel,X);
if ~exist(fullfile(Mopts.ExecutionPath,'myBeam.zip'),'file')
    pass = false ;
end
delete(fullfile(Mopts.ExecutionPath,'myBeam.mat')) ;
delete(fullfile(Mopts.ExecutionPath,'myBeam.zip')) ;

%% CASE 2: ARCHIVING IS SET TO 'IGNORE'. CHECK THAT THE TWO INPUT FILES AND THE OUTPUT FILE EXIST IN THE EXECUTION PATH
clear Mopts ;
Mopts.Type = 'UQLink' ;
Mopts.Name = 'my Beam 2' ;
Mopts.Command = [fullfile(uq_rootPath,'modules','uq_model','builtin','uq_uqlink','test','Case2', EXECNAME), ' SSB_Input_v2_1.inp', ' SSB_Input_v2_2.inp'] ;    
Mopts.ExecutionPath = fullfile(uq_rootPath,'modules','uq_model','builtin','uq_uqlink','test','Case2') ;
Mopts.Template = {'SSB_Input_v2_1.inp.tpl', 'SSB_Input_v2_2.inp.tpl'} ;
Mopts.Output.Parser= 'uq_read_SSBeamDeflection' ;
Mopts.Output.FileName = 'output.out' ;
Mopts.Archiving.Action = 'ignore' ;
Mopts.Display = 'quiet' ;

myModel = uq_createModel(Mopts) ;
Yval = uq_evalModel(myModel,X);
if ~exist(fullfile(Mopts.ExecutionPath,'SSB_Input_v2_1000001.inp'),'file') || ...
        ~exist(fullfile(Mopts.ExecutionPath,'SSB_Input_v2_2000001.inp'),'file') || ...
        ~exist(fullfile(Mopts.ExecutionPath,'output.out'),'file')
    pass = false ;
end
delete(fullfile(Mopts.ExecutionPath,'myBeam2.mat')) ;
delete(fullfile(Mopts.ExecutionPath,'SSB_Input_v2_1000001.inp')) ;
delete(fullfile(Mopts.ExecutionPath,'SSB_Input_v2_2000001.inp')) ;
delete(fullfile(Mopts.ExecutionPath,'output.out')) ;


end
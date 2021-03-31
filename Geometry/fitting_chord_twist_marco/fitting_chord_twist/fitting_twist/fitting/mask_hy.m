close all
clear all
clc
format long

parpool('local',2)

ObjectiveFunction = @mask_run;
nvars = 6;

n = 20;

B = [6 6 6 6 6 6];

LB = B-n;
UB = B+n;

hybridopts = psoptimset('Display','iter','InitialMeshSize',1,'MaxFunEvals',Inf,'MaxIter',Inf);
hybridopts = psoptimset(hybridopts,'UseParallel','always','CompletePoll','on','Vectorized','off');
hybridopts = psoptimset(hybridopts,'OutputFcns',@hyoutfun);

options = gaoptimset('Display','iter','Generations',Inf,'PopulationSize',220,'StallGenLimit',15,'TolFun',0);
options = gaoptimset(options,'UseParallel','always','Vectorized','off');
options = gaoptimset(options,'OutputFcns',@outfun);
options = gaoptimset(options,'HybridFcn',{@patternsearch,hybridopts});

dos('rmdir /S /Q .\History');
dos('rmdir /S /Q .\hyHistory');
dos('rmdir /S /Q .\Results');

dos('mkdir History');
dos('mkdir hyHistory');
dos('mkdir Results');

[x,fval,exitflag,output] = ga(ObjectiveFunction,nvars,[],[],[],[],LB,UB,[],options);
% [x,fval,exitflag,output] = patternsearch(ObjectiveFunction,B',[],[],[],[],LB,UB,[],hybridopts);

save ('.\Results\x.dat', 'x', '-ascii', '-double');
save ('.\Results\fval.dat', 'fval', '-ascii', '-double');
save ('.\Results\exitflag.dat', 'exitflag', '-ascii', '-double');

delete(gcp)
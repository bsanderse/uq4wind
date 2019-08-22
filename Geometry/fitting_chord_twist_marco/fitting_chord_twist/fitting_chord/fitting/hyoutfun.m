function [stop,options,optchanged] = hyoutfun(optimvalues,options,~)

stop = false;
optchanged = false;

save ('./hyHistory/x.dat', '-struct', 'optimvalues', 'x', '-ascii', '-double','-append');
save ('./hyHistory/fval.dat', '-struct', 'optimvalues', 'fval', '-ascii', '-double','-append');


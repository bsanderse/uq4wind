function [state,options,optchanged] = outfun(options,state,~)

optchanged = false;

save (['./History/population_generation_',int2str(state.Generation),'.dat'], '-struct', 'state', 'Population', '-ascii', '-double');
save (['./History/score_generation_',int2str(state.Generation),'.dat'], '-struct', 'state', 'Score', '-ascii', '-double');

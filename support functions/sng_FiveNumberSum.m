function [fns] = sng_FiveNumberSum(array)
%computes the five number summary, same as boxplot

Q = quantile(array,[.25 .5 .75]);
LF = Q(1) - 1.5 * (Q(3)-Q(1));
HF = Q(3) + 1.5 * (Q(3)-Q(1));

fns = [LF,Q,HF];

end
function [df] = df1dim(x, F, stimulus, P)
% Must accompany minfunc.

global pcom; % Defned in minfunc.
global xicom;
global nrdfun;

df  = dot(feval(nrdfun, pcom + x.*xicom, F, stimulus, P),xicom);

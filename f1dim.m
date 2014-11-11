% Move pcom (x units in xicom direction), and then evaluate the function there.

function [f] = f1dim(x, F, stimulus, P)
% Must accompany minfunc.

global pcom; % Defned in minfunc.
global xicom;
global nrfunc;

f  = feval(nrfunc, pcom + x.*xicom, F, stimulus, P);





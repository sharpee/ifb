% Move pcom (x units in xicom direction), and then evaluate the function there.

function [f] = f1dim(x, F, stimulus, P)
% Based on f1dim.c from Numerical Recipes in C (1992)
% by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery 

global pcom; % Defined in minfunc.
global xicom;
global nrfunc;

f  = feval(nrfunc, pcom + x.*xicom, F, stimulus, P);





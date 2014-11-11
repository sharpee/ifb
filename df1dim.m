function [df] = df1dim(x, F, stimulus, P)
% Based on df1dim.c from Numerical Recipes in C (1992)
% by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery 

global pcom; % Defined in minfunc.
global xicom;
global nrdfun;

df  = dot(feval(nrdfun, pcom + x.*xicom, F, stimulus, P),xicom);

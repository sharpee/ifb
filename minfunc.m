function [p, xi, fret] = minfunc(p, xi, func, dfunc, F, stimulus, P)

TOL    = 2.0e-5;
global pcom xicom nrfunc nrdfun;
nrfunc = func;
nrdfun = dfunc;
pcom   = p;
xicom  = xi;
ax     = 0.0;
xx     = 1;
ftemp=@f1dim;
ftemp2=@df1dim;
[ax, xx, bx, fa1, fc1, fb1] = bracket(ax, xx, ftemp, F, stimulus, P);
[fret, xmin] = minfinder(ax,xx,bx,ftemp,ftemp2,TOL, F, stimulus, P);
xi     = xi.*xmin;
p      = p + xi;

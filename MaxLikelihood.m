function [pfinal,fret] = MaxLikelihood(p, func, dfunc, F, ftol, stimulus, P)
% Based on frprmn.c from Numerical Recipes in C (1992)
% by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery 

ITMAX = 1000;
EPS   = 1.0e-10;
fp = feval(func, p, F, stimulus, P);
xi = feval(dfunc, p, F, stimulus, P);
g  = -xi;
h  = g;
xi = g;
fret_vec = zeros(ITMAX+1,1);
fret_vec(1) = fp;
figure; plot(fret_vec(1),'r-'); drawnow;
for its=1:ITMAX,
    disp(['Iteration ' num2str(its)]);
    [p, xi, fret] = minfunc(p, xi, func, dfunc, F, stimulus, P);
    fret_vec(its+1) = fret;
    plot(fret_vec(1:its+1),'r-'); drawnow;
    if  ( 2.0*abs(fret-fp) < ftol*( abs(fret) + abs(fp) + EPS ))
        pfinal = p;
        break;
    end
    if sum(isnan(p))>0
        pause(1)
    end
    fp = fret;
    xi = feval(dfunc, p, F, stimulus, P);
    gg = sum(g.^2);
    dgg = sum( (xi + g).*xi );
    if gg == 0,
        break;
    end
    gam = dgg/gg;
    g = -xi;
    h = g + gam.*h;
    xi = h;
end

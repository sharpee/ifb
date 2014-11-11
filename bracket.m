function [ax, bx, cx, fa, fb, fc]=bracket(ax,bx,f, F, stimulus, P)
% Based on bracket.c from Numerical Recipes in C (1992)
% by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery 

GOLD=1.618034;
GLIMIT=100.0;
TINY=1.0e-20;
fa=f(ax, F, stimulus, P);
fb=f(bx, F, stimulus, P);
if fb>fa,
    dum=ax ; ax=bx ; bx=dum ;
    dum=fb ; fb=fa ; fa=dum ;
end;
cx=bx+GOLD*(bx-ax);
fc=f(cx, F, stimulus, P);
while fb>fc
    r=(bx-ax)*(fb-fc);
    q=(bx-cx)*(fb-fa);
    u= bx -((bx-cx)*q-(bx-ax)*r)/(2.0*max([abs(q-r),TINY])*sign(q-r));
    ulim=(bx)+GLIMIT*(cx-bx);
    if (bx-u)*(u-cx)>0.0,
        fu=f(u, F, stimulus, P) ;
        if (fu < fc)
            ax=bx;
            bx=u;
            fa=fb;
            fb=fu;
            return;
        elseif (fu> fb)
            cx=u;
            fc=fu;
            return;
        end
        u=cx+GOLD*(cx-bx);
        fu=f(u, F, stimulus, P);
    elseif (cx-u)*(u-ulim)>0.0,
        fu=f(u, F, stimulus, P) ;
        if fu<fc
            bx=cx;
            cx=u;
            %u=u+GOLD*(u-cx);
            u = cx+GOLD*(cx-bx);
            fb=fc;
            fc=fu;
            fu=f(u, F, stimulus, P);
        end
    elseif ((u-ulim)*(ulim-cx))>=0.0
        u=ulim;
        fu=f(u, F, stimulus, P);
    else
        u=(cx)+GOLD*(cx-bx);
        fu=f(u, F, stimulus, P) ;
    end ;
    ax=bx ; bx=cx ; cx=u ;
    fa=fb ; fb=fc ; fc=fu ;
end




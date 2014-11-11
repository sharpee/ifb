function [rval, xmin] = minfinder(ax, bx, cx, f, df, tol, F, stimulus, P)

ITMAX = 500;
ZEPS  = 1.0e-9;
e     = 0.0;
if ax < cx,
   a = ax;
   b = cx;
else
   a = cx;
   b = ax;
end
x  = bx;
w  = bx;
v  = bx;
fx = feval(f,x, F, stimulus, P);  fw = fx; fv = fx;
dx = feval(df,x, F, stimulus, P); dw = dx; dv = dx;
for iter=1:ITMAX,
   xm   = 0.5*(a+b);
   tol1 = tol*abs(x)+ZEPS;
   tol2 = 2.0*tol1;
   if abs(x-xm) <= (tol2-0.5*(b-a))
      xmin  = x;
      rval  = fx;
      return;
   end
   if (abs(e) > tol1)
      d1 = 2.0*(b-a);
      d2 = d1;
      if (dw ~= dx),
         d1 = (w-x)*dx/(dx-dw);
      end
      if (dv ~= dx),
         d2 = (v-x)*dx/(dx-dv);
      end
      u1   = x+d1;
      u2   = x+d2;
      ok1  = ((a-u1)*(u1-b) > 0.0) & (dx*d1 <= 0.0);
      ok2  = ((a-u2)*(u2-b) > 0.0) & (dx*d2 <= 0.0);
      olde = e;
      e    = d;
      if ok1 | ok2,
         if ok1 & ok2,
            if abs(d1) < abs(d2), d = d1; else d = d2; end
         elseif ok1,
            d = d1;
         else
            d = d2;
         end;

         if abs(d) <= abs(0.5*olde),
            u=x+d;
            if (u-a < tol2) | (b-u < tol2),
               d = tol1*signNR(xm-x);
            end
         else
            if dx >= 0, e=a-x; else e=b-x; end
            d = 0.5*e;
         end
      else
         if dx >= 0, e=a-x; else e=b-x; end
         d = 0.5*e;
      end
   else
      if dx >= 0, e=a-x; else e=b-x; end
      d = 0.5*e;
   end
   if abs(d) >= tol1,
      u  = x + d;
      fu = feval(f,u, F, stimulus, P);
   else
      u  = x + tol1*signNR(d);
      fu = feval(f,u, F, stimulus, P);
      if fu > fx,
         xmin = x;
         rval = fx;
         return;
      end
   end
   du = feval(df, u, F, stimulus, P);
   if fu <= fx,
      if u >= x, a=x; else b=x; end;
      v = w; fv = fw; dv = dw;
      w = x; fw = fx; dw = dx;
      x = u; fx = fu; dx = du;
   else
      if u < x, a=u; else b=u; end;
      if (fu <= fw) | (w == x),
         v = w; fv = fw; dv = dw;
         w = u; fw = fu; dw = du;
      elseif (fu < fv | v == x | v == w),
         v = u; fv = fu; dv = du;
      end;
   end;
end;
xmin = x;
rval = fx;

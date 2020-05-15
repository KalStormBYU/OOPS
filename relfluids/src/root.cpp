#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include <root.h>

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define LTRACE 0

/*-------------------------------------------------------------------------
 *
 *  func - returns f(x0), df/dx(x0).  The arguments are 
 *         func(f, df, x0, fpar)
 *
 *-------------------------------------------------------------------------*/
int rtsafe(int (*func)(double *, double *, double , double *), 
              double *rts, double x1, double x2, double tol, 
              double *fpar,int trace)
{

  int j, rc;
  double df, dx, dxold, f, fh, fl;
  double tmp, xh, xl;

  const double MAXIT  = 50;

  rc = (*func)(&fl, &df, x1, fpar);
  if (rc != 1) return -1;
  rc = (*func)(&fh, &df, x2, fpar);
  if (rc != 1) return -1;

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    if ( trace ) {
      printf("RTSAFE: Error.  Root must be bracketed.\n");
      printf("x1 = %g, f1 = %g\n",x1,fl);
      printf("x2 = %g, f2 = %g\n",x2,fh);
    }
    return -1;
  }

  if (fl == 0.0) {
    *rts = x1;
    return 1;
  }
  if (fh == 0.0) {
    *rts = x2;
    return 1;
  }

  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  }
  else {
    xl = x2;
    xh = x1;
  }

  *rts = 0.5*(x1 + x2);
  dxold = fabs(x2 - x1);
  dx = dxold;
  rc = (*func)(&f, &df, *rts, fpar);
  if (trace) {
    printf("RTSAFE: rts=%g, f=%g, df=%g\n",*rts, f, df);
  }

  if (rc != 1) return -1;

  for (j=0; j < MAXIT; j++) {
    if ((((*rts-xh)*df-f)*((*rts-xl)*df-f) > 0.0) || 
          (fabs(2.0*f) > fabs(dxold*df))) {
      /* Bisect if Newton solve would give root out of range, or if the
       * solver isn't converging quickly.
       */
      dxold = dx;
      dx = 0.5*(xh-xl);
      *rts = xl + dx;
      if (xl == *rts) return 1;
    }
    else {
      /* Newton solve */
      dxold = dx;
      dx = f/df;
      tmp = *rts;
      *rts -= dx;
      if (tmp == *rts) return 1;
    }
    if (fabs(dx) < tol) return 1;

    rc = (*func)(&f, &df, *rts, fpar);
    if (trace) {
      printf("RTSAFE: rts=%g, f=%g, df=%g\n",*rts, f, df);
    }
    if (rc != 1) return -1;

    if (f < 0.0)
      xl = *rts;
    else
      xh = *rts;
  }

  if ( trace ) {
    printf("RTSAFE:  No solution.\n");
  }
  return -1;

}


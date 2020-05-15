#ifndef ROOT_H
#define ROOT_H

int rtsafe(int (*func)(double *, double *, double, double *),
              double *rts, double x1, double x2, double tol,
              double *fpar, int trace);

#endif

#include <relfluid.h>
#include <operators.h>
#include <root.h>
#include <iostream>
#include <cmath>

//  Constructor
Fluid::Fluid(Domain &d, Solver &s) : ODE(6,0)
{
    if (d.getGhostPoins() < 2)
    {
std:cerr << "Warning: Domain has fewer ghost points than expected" << endl;
    }
    domain = &d;
    solver = &s;

    // Set parameters
    params = new Parameters();

    reallocateData();
}

// Destructor
Fluid::~Fluid()
{
    delete params;
}


// Initial data

Fluid::initData()
{
    double x0 = 0.0;
    double sigma = 0.5;

    for (auto it = data.begin(); it != data.end(); it++)
    {
        const double *x = it->getGrid().getPoint();
        unsigned int nx = it->getGrid().getSize();
        double **u = it->getData();
        for (unsigned int i = 0; i < nx; i++)
        {
            u[V_RHO][i] = std::exp(-(x[i]-x0)*(x[i]-x0)/(sigma*sigma)) + floor;
            vv[V_VX][i] - 0.0;
            v[V_P][i] = k * pow(v[V_RHO][i], Gamma);
        }
    }
}


// SOURCE
// primToCon {{{
void Fluid::primToCon(){
  for(auto it = data.begin(); it != data.end(); ++it){
    const double *x = it->getGrid().getPoints();
    unsigned int nx = it->getGrid().getSize();
    double **u = it->getData();
    for(unsigned int i = 0; i < nx; i++){
      double W = calcW(u[V_VX][i]);
      double h = calcH(u[V_RHO][i], u[V_P][i]);

      u[U_D  ][i] = u[V_RHO][i]*W;
      u[U_SX ][i] = h*W*W*u[V_VX][i];
      u[U_TAU][i] = h*W*W - u[V_P][i] - u[U_D][i];
    }
  }
}
// }}}

// calcW {{{
double Fluid::calcW(double v){
  return 1.0/std::sqrt(1.0 - v*v);
}
// }}}

// calcH {{{
double Fluid::calcH(double rho, double P){
  FluidParameters *fp = (FluidParameters*) params;
  double gamma = fp->getGamma();

  return rho + gamma/(gamma - 1.0)*P;
}
// }}}

// doAfterBoundaries {{{
void Fluid::doAfterBoundaries(bool intermediate){
  // Convert the conserved variables to primitive variables.
  unsigned int nb = domain->getGhostPoints();
  for(auto it = data.begin(); it != data.end(); ++it){
    const double *x = it->getGrid().getPoints();
    unsigned int nx = it->getGrid().getSize();
    double **u;
    if(!intermediate){
      u = it->getData();
    }
    else{
      u = it->getIntermediateData();
    }
    
    double upt[nEqs];
    for(unsigned int i = nb; i < nx - nb; i++){
      // Grab a single data point for the primitive solver.
      for(unsigned int m = 0; m < nEqs; m++){
        upt[m] = u[m][i];
      }
      // Run the primitive solver.
      conToPrimPt(upt);

      // Return the point data.
      for(unsigned int m = 0; m < nEqs; m++){
        u[m][i] = upt[m];
      }
    }
  }
}
// }}}

// conToPrimPt {{{
void Fluid::conToPrimPt(double *upt){
  // Get some frequently used parameters.
  // Ryan, TJ, and Elijah: assuming your floor, gamma, and so forth
  // are defined in your fluid class, you can get rid of all this stuff.
  // If you named your floor, gamma, and so forth differently than I did,
  // you'll need to go through this function to update them.
  FluidParameters *fp = (FluidParameters*) params;
  double floor = fp->getFloor();
  double gamma = fp->getGamma();


  // Apply the floor.
  // FIXME  The flooring given below depends on D and tau being below
  // vacuum values independently.  This should probably be imposed in
  // a more intelligent way that considers more of the possibilities.
  if(upt[U_D] < floor){
    upt[U_D] = floor;
    upt[U_SX] = 0.0;
    upt[U_TAU] = floor;
  }
  else{
    if (upt[U_TAU] < floor){
      upt[U_TAU] = floor;
    }

    double Ssq = upt[U_SX]*upt[U_SX];
    double Ssq_max = (2.0*upt[U_D] + upt[U_TAU])*upt[U_TAU];
    if(Ssq_max < Ssq){
      double t_min = VELOCITY_FACTOR * std::sqrt(Ssq_max/Ssq);
      upt[U_SX] *= t_min;
    }
  }

  // Call the EOS solver.
  // Set up some work variables for the EOS solver.
  double fpar[NFPAR];
  conToPrimWorkVars(fpar, upt, gamma, floor);

  int rc = MHDSolver4(upt, fpar);
  
  // Check that pressure is still physical.
  if(upt[V_P] < floor){
    upt[V_P] = floor;
  }
}
// }}}

// conToPrimWorkVars {{{
void Fluid::conToPrimWorkVars(double *fpar, double *u, double gamma, double vacuum){
  double kappa;

  double vsq = u[V_VX]*u[V_VX];

  if(u[V_VX]*u[V_VX] >= 1.0){
    std::cout << "  conToPrimWorkVars:  The velocity is superluminal.  To correct \n";
    std::cout << "  conToPrimWorkVars:  this, we simply set it to 0.99999999 \n";
    std::cout << "  conToPrimWorkVars:  which seems sort of stupid. \n";
    std::cout << "  conToPrimWorkVars:       vsq = " << (u[V_VX]*u[V_VX]) << " \n";
    //FIXME It would be well to come up with something more intelligent
    vsq = 0.99999999;
  }

  // Define some temporary and auxiliary variables.
  fpar[FP_SSQ]    = u[U_SX]*u[U_SX];
  fpar[FP_D]      = u[U_D];
  fpar[FP_TAU]    = u[U_TAU];
  fpar[FP_GAMMA]  = gamma;
  fpar[FP_VSQ]    = vsq;
  fpar[FP_COORD_X] = 0.0; // FIXME: We should pass in the x coordinate.
  kappa            = u[V_P]/std::pow(u[V_RHO],gamma);
  fpar[FP_KAPPA]   = ( kappa < 1000.0 ) ? kappa : 1000.0;

  fpar[FP_TRACE]   = 0.0;

  fpar[FP_C2PWARN] = 1;
  fpar[FP_VACUUM]  = vacuum;
}
// }}}

// MHDSolver4 {{{
int Fluid::MHDSolver4(double *u, double *fpar){
  double gamma = fpar[FP_GAMMA];
  double floor = fpar[FP_VACUUM];

  double Ssq = fpar[FP_SSQ];
  double D   = u[U_D];
  double tau = u[U_TAU];
  double Sx  = u[U_SX];

  int j;
  double qf_1, qf_2, qg_2, qf_max;
  double ql, qh, f_at_qh, f_at_qh_v2;
  double f_at_ql = 0.0;
  double rescale_factor1;
  double a = tau + D;
  double t1 = 0.5 * gamma;
  double beta_sq = Ssq/(a*a);
  double t2 = t1*t1 - (gamma-1.0)*beta_sq;
  // t2 shows up under a square root, so we check to see if it is < 0.
  if (t2 < 0.0){
    std::cout << "MHDSolver4:  problem with upper bound. t2 < 0. t2 = " << t2 << "\n";
    exit(2);
  }
  double sqrt_t2 = std::sqrt(t2);
  qf_1 = t1 - sqrt_t2; // The smaller root of f
  qf_2 = t1 + sqrt_t2; // The larger root of f
  qf_max = t1;         // The location of the max of f
  qg_2 = std::sqrt(beta_sq); // The larger root of g

  if (std::isnan(qg_2)){
    std::cout << "qg_2 is a nan.\n";
  }

  /*
  If the velocity vector (or S^i) is the zero vector, the conserved to
  primitive inversion is analytic.  This translates to beta_sq being zero.
  We do this case separately and return.
  */
  if (beta_sq < 1.e-16){
    u[V_RHO] = D;
    u[V_P  ] = (gamma - 1.0)*tau;
    u[V_VX ] = 0.0;

    return 1;
  }

  if ( fabs( (qg_2 - qf_1) / (qg_2 + qf_1 + 1.0e-15) ) < 1.e-6 ){
    qg_2 += 1.e-5;
  }

  // Set the lower bound of qg_2 unless qg_2 is smaller than the location of
  // the maximum of f.
  ql = fmax( qg_2, qf_max );

  //  This temporarily initializes the upper bound to a possible value.
  //  The following code checks this to make sure it works.
  qh = qf_2 ;
  /* {{{ ... documentation on the cases that need to be considered ...                                                 We consider various possibilities for the relative location of qg_2 and
  qf_2 as they should be the lower and upper bounds, respectively, for our
  root, q, of F(q).

  Case 1 is that qg_2 < qf_2 and a physical inequality is obeyed.  All is well
  in this case and the Newton solve should proceed without issue.
  
  Case 2 considers the violation of another physical inequality and adjusts
  (rescales) the components of S in order to enforce the inequality.
  
  Case 3 considers the numerical "vanishing" of D/a relative to beta_sq which
  can also numerically violate physical inequalities.

  Case 4 considers the (seemingly remote but nonetheless realized) possibility
  that qg_2 and qf_2 are exactly the same (to machine precision).

  No claim is made that these cases are exhaustive.  They were simply ones
  that arose in debugging the code and trying to make this primitive solver as
  robust as possible.  There is a sort of "ad hoc-ness" to some of these.
  Could they be improved?  Probably.  But for now they seem to work.
  }}}  */
  // CASE 1
  if ( qg_2 < qf_2  &&  beta_sq < 1.0 ){
    qh = qf_2;
  }
  // CASE 2
  else if ( beta_sq + (D/a)*(D/a) - 1.0 >= 0.0 ){
    /*
    We rescale the components of Sd (and hence Ssq) in order that
    the inequality is enforced.

    FIXME:  Note that this scaling factor is a bit arbitrary.  It would
    FIXME:  be well to experiment with the value added in.
    */
    rescale_factor1 = 1.0 / (std::sqrt(beta_sq + (D/a)*(D/a)) + 1.e-14);

    // Now multiply Sx by the scaling factor.
    Sx *= rescale_factor1;

    // Recalculate Ssq.
    Ssq = Sx*Sx;

    // Now we have to recalculate the roots of f and g with this new Ssq.
    beta_sq = Ssq/(a*a);
    t2 = t1*t1 - (gamma-1.0)*beta_sq;
    if (t2 < 0.0){
      std::cout << "  MHDSolver4: t2 < 0, i.e. t2=" << t2 << "\n";
    }
    qf_2 = t1 + std::sqrt(t2); // new root of f; upper bound on q
    qg_2 = std::sqrt(beta_sq); // new root of g; lower bound on q

    // Check that these new roots are inded good upper and lower bounds.
    if ( qg_2 < qf_2 ){
      ql = qg_2;
      qh = qf_2;
    }
    else{
      
    }

    // Because we made a rescaling of the components of S and of Ssq, we need
    // to save these.
    u[U_SX] = Sx;
    fpar[FP_SSQ] = Ssq;
  }
  // CASE 3
  else if ( fabs( (beta_sq + (D/a)*(D/a)) - beta_sq) <= 6.0e-16 ){
    // {{{
    std::cout << "  MHDSolver4:  SURPRISE! I deleted the treatment of this \n";
    std::cout << "  MHDSolver4:  possibility as I didn't think it ever arose. \n";
    std::cout << "  MHDSolver4:  I guess it does after all. \n";
    // }}}
  }
  // CASE 4
  else if (qg_2 == qf_2){
    qg_2 -= 1.e-16;
    qf_2 += 1.e-16;

    ql = qg_2;
    qh = qf_2;
  }
  else{
    //{{{
    // FIXME:
    // This else is effectively a garbage dump for all the other possibilities
    // not delineated in the if/else above.  It is a sort of "none of the 
    // above".  As a result, if we get here, we essentially don't know why ...
    std::cout << "  MHDSolver4:    something unexpected (1) \n"
              << "       qg_2 = " << qg_2 << "\n"
              << "       qf_2 = " << qf_2 << "\n"
              << "    beta_sq = " << beta_sq << "\n"
              << "       D    = " << D << "\n"
              << "       tau  = " << tau << "\n"
              << "       Ssq  = " << Ssq << "\n"
              << "    beta_sq + (D/a)^2 = " << beta_sq+(D/a)*(D/a) << "\n"
              << "   dump fpar  \n"
              << "      fpar[FP_SSQ]     = " << fpar[FP_SSQ] << "\n"
              << "      fpar[FP_D  ]     = " << fpar[FP_D  ] << "\n"
              << "      fpar[FP_TAU]     = " << fpar[FP_TAU] << "\n"
              << "      fpar[FP_GAMMA]   = " << fpar[FP_GAMMA] << "\n"
              << "      fpar[FP_VSQ]     = " << fpar[FP_VSQ] << "\n"
              << "      fpar[FP_COORD_X] = " << fpar[FP_COORD_X] << "\n";
    //}}}
  }


  // At this point, we should have the two values, namely qg_2 and qf_2,
  // that should form a legitimate bracket for finding the root of F(q).
  // Let's check.
  //
  // We want to check that at ql, the function F(q) := f(q)-(const)*sqrt(g(q))
  // is positive.  First just calculate F(ql).
  int rc = 0;
  rc = func4_p(&f_at_ql, ql, fpar);
  if (rc < 0){
    std::cout << "  MHDSolver4:  Something's amiss with ql \n";
    std::cout << "  MHDSolver4:  We couldn't calculate f(ql). \n";
  }

  // Now check that at qh, the function F(q) := f(q)-(const)*sqrt(g(q))
  // is negative.  For now, just calculate F(qh).
  rc = 0;
  rc = func4_p(&f_at_qh, qh, fpar);
  if (rc < 0){
    std::cout << "  MHDSolver4:  Something's amiss with qh \n";
    std::cout << "  MHDSolver4:  We couldn't calculate f(qh). \n";
  }

  // If, indeed, F(ql) > 0 and F(qh) < 0, then we should have a
  // bracket for the root. We will check for the failure of this and do
  // what we can.
  if (f_at_ql < 0.0 || f_at_qh > 0.0){
  
  }

  if(f_at_ql > 0.0 && f_at_qh > 0.0){
    qh += 6.e-16;

    rc = func4_p (&f_at_qh_v2, qh, fpar) ;
  }

  /*
  At this point, we should have a legitimate bracket for our root, namely
  [ql, qh] := [qg_2, qf_2].  We should now be able to proceed with the 
  Newton solve.

  As our initial guess for the root, take the midpoint of the bracket.
  */
  double qi = 0.5*(ql + qh);
  double q = qi;

  //  Incorporate a little reduncancy in the event things go bad.
  double ql0 = ql;
  double qh0 = qh;
  double q0 = q;

  //  Set trace parameter for rtsafe to 0 (no tracing), set the tolerance for
  //  the Newton solve, and call rtsafe to get the root, q.
  int xtrace = 0;
  double xtol = 1.0e-14;
  int xrc = rtsafe(func4_p_d, &q, ql, qh, xtol, fpar, xtrace);

  // In the event of failure, a postmortem ...
  if (xrc < 0){
    std::cout << "  MHDSolver4:  The Newton solve (rtsafe) failed.\n"
              << "  MHDSolver4:  Current values include \n"
              << "    ql = " << ql << " \n"
              << "    qh = " << qh << " \n"
              << "    qi = " << qi << " \n"
              << "    q   = " << q << " \n"
              << "    ql0 = " << ql0 << " \n"
              << "    qh0 = " << qh0 << " \n"
              << "    q0  = " << q0 << " \n"
              << "  MHDSolver4:  Some conserved quantities and \n"
              << "  MHDSolver4:  inequalities, \n"
              << "    D               = " << D << " \n"
              << "    tau+D           = " << a << " \n"
              << "    D^2/a^2         = " << (D*D)/(a*a) << " \n"
              << "    Ssq/(a*a)       = " << Ssq/(a*a) << " \n"
              << "    (Ssq+D^2)/(a*a) = " << (Ssq+D*D)/(a*a) << " \n"
              << "  MHDSolver4:  After rtsafe failure ... redoing ... \n";
    // Retrying rtsafe but with a different trace parameter so as to get
    // some more info out of rtsafe.
    ql = ql0;
    qh = qh0;
    q = q0;

    xtrace = 1;
    xrc = rtsafe(func4_p_d, &q, ql, qh, xtol, fpar, xtrace);
    exit(2);
    return -4;
  }

  //  If we get to this point we may actually have a solution.
  //  We have to rescale back by a = tau+D;
  q *= a;

  // Knowing the value for q (:= hW^2, i.e. related to the enthalpy), we
  // can now get the remaining primitive variables.  Start with the velocity
  // squared: v^2.
  double vsq = Ssq / (q*q);

  // Now calculate the primitive variables.
  u[V_RHO] = D*std::sqrt(1.0 - vsq);
  u[V_P  ] = (gamma-1.0)/gamma*(q*(1.0-vsq) - u[V_RHO]);
  u[V_VX ] = Sx / q;

  if(u[V_P] == 0){
    u[V_P] = (gamma - 1.0)*tau;
  }

  if(u[V_RHO] < floor){
    u[V_RHO] = floor;
    u[V_VX] = 0.0;
  }

  if (u[V_P] < floor){
    u[V_P] = floor;
  }
  assert( D > 0.0 );

  return 0;
}
// }}}

// func4_p {{{
int Fluid::func4_p(double *f, double q, double *fpar){
  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double tau   = fpar[FP_TAU];
  double gamma = fpar[FP_GAMMA];

  double dis;

  double a = tau+D;
  dis = q*q - Ssq/(a*a);
  if ( dis < 1.e-18){
    dis = 1.e-18;
  }

  *f  = -(q - 0.5*gamma)*(q - 0.5*gamma) + gamma*gamma/4.0
      - (gamma-1.0)*Ssq/(a*a) - (gamma-1.0)*D/a*std::sqrt(dis);

  return 1;
}
// }}}

// func4_p_d {{{
int Fluid::func4_p_d(double *f, double *df, double q, double *fpar){
  double Ssq   = fpar[FP_SSQ];
  double D     = fpar[FP_D];
  double tau   = fpar[FP_TAU];
  double gamma = fpar[FP_GAMMA];

  double dis;

  double a = tau+D;
  dis = q*q - Ssq/(a*a);
  if (dis < 1.e-18){
    dis = 1.e-18;
  }

  *f  = - (q - 0.5*gamma)*(q - 0.5*gamma) + gamma*gamma/4.0
        - (gamma-1.0)*Ssq/(a*a) - (gamma-1.0)*D/a*std::sqrt(dis);

  *df = -2.0*(q - 0.5*gamma) - (gamma-1.0)*D/a*q*std::sqrt(dis);

  return 1;
}
// }}}

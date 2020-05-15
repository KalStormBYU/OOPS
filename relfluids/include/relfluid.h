#ifndef RELFLUID_H
#define RELFLUID_H

#include <ode.h>

class Fluid : public ODE {
    private:
        // Variable labels
        static const unsigned int U_D   = 0;
        static const unsigned int U_SX  = 1;
        static const unsigned int U_TAU = 2;
        static const unsigned int V_RHO = 3;
        static const unsigned int V_VX  = 4;
        static const unsigned int V_P   = 5;
        const double floor              = 1.0e-6;
        const double k                  = 1.0;
        const double Gamma              = 5/3;
        // Conversions between primitive and conserved variables.
        void primToCon();
        void conToPrimPt(double *u);

        // Utility functions for the primitive solver
        void conToPrimWorkVars(double *fpar, double *u, double gamma, double vacuum);
        int MHDSolver4(double *upt, double *fpar);
        static int func4_p(double *f, double q, double *fpar);
        static int func4_p_d(double *f, double *df, double q, double *fpar);

        // Some utility variables for the priitive solver
        static constexpr double VELOCITY_FACTOR = 0.9999;
        static const unsigned int FP_SSQ      = 0;
        static const unsigned int FP_D        = 1;
        static const unsigned int FP_TAU      = 2;
        static const unsigned int FP_GAMMA    = 3;
        static const unsigned int FP_VSQ      = 4;
        static const unsigned int FP_KAPPA    = 5;
        static const unsigned int FP_TRACE    = 6;
        static const unsigned int FP_C2PWARN  = 7;
        static const unsigned int FP_VACUUM   = 8;
        static const unsigned int FP_COORD_X  = 9;
        static const unsigned int NFPAR       = 10;

        // Useful utility functions
        double calcW(double v);
        double calcH(double rho, double P);

    protected:
        virtual void applyBoundaries (bool intermediate);
        virtual void rhs (const Grid& grid, double **u, double **dudt);

    public:
        Fluid(Domain &d, Solver &s);
        virtual ~Fluid();
        virtual void initData();
        virtual void doAfterBoundaries(bool intermediate);

}


#endif

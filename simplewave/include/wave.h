#ifndef WAVE_H
#define WAVE_H

#include <ode.h>

class Wave : public ODE
{
    private:
        // Variable labels
        static const unsigned int U_PHI  = 0;
        static const unsigned int U_PI   = 1;
        static const unsigned int U_CHI  = 2;

    protected:
        virtual void applyBoundaries(bool intermediate);
        virtual void rhs(const Grid& grid, double **u, double **dudt);

    public:
        Wave(Domain& d, Solver& s);
        virtual ~Wave();

        virtual void initData();
};

#endif

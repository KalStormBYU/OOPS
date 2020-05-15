#include <wave.h>
#include <operators.h>
#include <iostream>
#include <waveparameters.h>
#include <cmath>

// Constructor
Wave::Wave(Domain& d, Solver& s) : ODE(2,1)
{
    if (d.getGhostPoints() < 2) 
    {
        std::cerr << "Warning: domain has fewer ghost points than expected. Expect incorrrect behavior. \n";
    }
    domain = &d;
    solver = &s;

    // Set some default parameters
    params = new Parameters();
    reallocateData();
}

// Destructor
Wave::~Wave()
{
    // delete params;
}


void Wave::rhs(const Grid& grid, double **u, double **dudt)
{
    // Define things we need
    double stencil3[3] = {0.0 , 0.0 , 0.0};
    double stencil5[5] = {0.0 , 0.0 , 0.0 , 0.0 , 0.0};
    double dx = grid.getSpacing();
    int shp = grid.getSize();
    // Calculate the left boundary. We switch to a different operator on the boundaries, which should
    // just be ghost points that will be overwritten anyway.
    // Leftmost Point
    dudt[U_PHI][0] = u[U_PI][0];
    stencil3[0] = u[U_CHI][0];
    stencil3[1] = u[U_CHI][1];
    stencil3[2] = u[U_CHI][2];
    dudt[U_PI][0] = operators::dx_2off(stencil3, dx);
    stencil3[0] = u[U_PI][0];
    stencil3[1] = u[U_PI][1];
    stencil3[2] = u[U_PI][2];
    dudt[U_CHI][0] = operators::dx_2off(stencil3, dx);

    // Second leftmost point
    dudt[U_PHI][1] = u[U_PI][1];
    stencil3[0] = u[U_CHI][0];
    stencil3[1] = u[U_CHI][1];
    stencil3[2] = u[U_CHI][2];
    dudt[U_PI][1] = operators::dx_2(stencil3, dx);
    stencil3[0] = u[U_PI][0];
    stencil3[1] = u[U_PI][1];
    stencil3[2] = u[U_PI][2];
    dudt[U_CHI][1] = operators::dx_2(stencil3, dx);

    // No all interior points
    for (int i = 2; i < shp - 2; i++)
    {
        dudt[U_PHI][i] = u[U_PI][i];

        for (int j = 0; j < 5; j++)
        {
            stencil5[j] = u[U_CHI][i-2+j];
        }
        dudt[U_PI][i] = operators::dx_4(stencil5, dx);

        for (int j = 0; j < 5; j++)
        {
            stencil5[j] = u[U_PI][i-2+j];
        }
        dudt[U_CHI][i] = operators::dx_4(stencil5,dx);
    }

    // Second rightmost point
    dudt[U_PHI][shp - 2] = u[U_PI][shp - 2];
    stencil3[0] = u[U_CHI][shp - 3];
    stencil3[1] = u[U_CHI][shp - 2];
    stencil3[2] = u[U_CHI][shp - 1];
    dudt[U_PI][shp - 2] = operators::dx_2(stencil3, dx);
    stencil3[0] = u[U_PI][shp - 3];
    stencil3[1] = u[U_PI][shp - 2];
    stencil3[2] = u[U_PI][shp - 1];
    dudt[U_CHI][shp - 2] = operators::dx_2(stencil3, dx);

    // Rightmost point
    dudt[U_PHI][shp - 1] = u[U_PI][shp - 1];
    stencil3[0] = u[U_CHI][shp - 3];
    stencil3[1] = u[U_CHI][shp - 2];
    stencil3[2] = u[U_CHI][shp - 1];
    dudt[U_PI][shp - 1] = operators::dx_2off(stencil3, dx);
    stencil3[0] = u[U_PI][shp - 3];
    stencil3[1] = u[U_PI][shp - 2];
    stencil3[2] = u[U_PI][shp - 1];
    dudt[U_CHI][shp - 1] = operators::dx_2off(stencil3, dx);
}

void Wave::applyBoundaries(bool intermediate)
{
    unsigned int nb = domain->getGhostPoints();
    // Grab the data at the leftmost grid and the rightmostgrid
    auto left_it = data.begin();
    auto right_it = --data.end();

    double **left;
    double **right;

    if (intermediate)
    {
        left = left_it->getData();
        right = right_it->getData();
    }
    else
    {
        left = left_it->getIntermediateData();
        right = right_it->getIntermediateData();
    }
    unsigned int nr = right_it->getGrid().getSize();

    // Apply Neumann boundary conditions
    for (unsigned int i = 0; i < nb; i++)
    {
        left[U_PHI][i] = left[U_PHI][nb];
        left[U_PI][i] = left[U_PI][nb];
        left[U_CHI][i] = 0.0;

        right[U_PHI][nr - 1 - i] = right[U_PHI][nr-nb-1];
        right[U_PI][nr - 1 - i] = right[U_PI][nr-nb-1];
        right[U_CHI][nr - 1 - i] = 0.0;
    }
}

void Wave::initData()
{
    // Center of our Gaussian
    double x0 = 0.5;

    // loop through every grid and start assigning points
    for (auto it = data.begin(); it != data.end(); ++it)
    {
        const double *x = it->getGrid().getPoints();
        unsigned int nx = it->getGrid().getSize();
        double **u = it->getData();
        for (unsigned int i = 0; i < nx; i++)
        {
            double val = std::exp(-(x[i] - x0) * (x[i] - x0) * 64.0);
            u[U_PHI][i] = val;
            u[U_PI][i] = 0.0;
            u[U_CHI][i] = -128.0*(x[i] - x0)*val;

        }
    }

}

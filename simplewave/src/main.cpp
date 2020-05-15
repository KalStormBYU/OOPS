#include <domain.h>
#include <grid.h>
#include <rk4.h>
#include <cmath>
#include <cstdio>
#include <wave.h>
#include <waveparser.h>
#include <polynomialinterpolator.h>
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 2) 
    {
        std::cerr << "Usage: " << argv[0] << " <parameter file>" << std::endl;
        exit(1);
    }
    WaveParameters pars;
    WaveParser parser;
    parser.updateParameters(argv[1], &pars);

    char *fnames[1];
    fnames[0] = new char[16];
    sprintf(fnames[0], "phi");



    // Construct our domain and a grid to fit on it
    Domain domain = Domain();
    double bounds[2];
    bounds[0] = pars.getxmin();
    bounds[1] = pars.getxmax();
    int N = pars.getnpoints();
    domain.setBounds(bounds);
    double cfl = pars.getcfl();
    double tmax = pars.gettmax();
    domain.setCFL(cfl);
    int output_frequency = pars.getoutput_frequency();

    std::cout << "Creating grid with " << N << " points and bounds [" << bounds[0] << ", " << bounds[1] << "}" << std::endl;

    domain.addGrid(bounds, N);
    
    // Set up out ODE system
    RK4 rk4 = RK4();
    PolynomialInterpolator interpolator = PolynomialInterpolator(4);

    Wave ode = Wave(domain, rk4);    // Change to Fluid for Fluid ODE
    ode.setInterpolator(&interpolator);
    ode.initData();

    double ti = 0.0;
    double dt = domain.getCFL() * (--domain.getGrids().end())->getSpacing();
    unsigned int M = (tmax - ti)/dt;
    ode.output_frame(fnames[0], 0.0, 0);
    int it = 0;
//    ode.outputVTKScalar(fnames[0], 0.0,it,0);

    double t = 0.0;
    for (unsigned int i = 1; i <= M; i++)
    {
        ode.evolveStep(dt);
        ti += dt;
        if (i % output_frequency == 0)
        {
            std::cout << "Step " << i << " time " << t << std::endl;
            ode.output_frame(fnames[0], t, 0);
            it++;
//            ode.outputVTKScalar(fnames[0], t, it, 0);
        }
    }
/*    ode.dump_csv("phi00000.csv", 0, 0);
    for (unsigned int i = 0; i < M; i++)
    {
        double t = (i + 1)*dt;
        ode.evolveStep(dt);

        char buffer[12];
        sprintf(buffer, "phi%05d.csv", i+1);
        ode.dump_csv(buffer, t, 0);
    }
    */
    delete [] fnames[0];
    return 0;
}

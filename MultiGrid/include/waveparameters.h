#ifndef WAVE_PARAMETERS_H
#define WAVE_PARAMETERS_H

#include <parameters.h>

#include <string>

// DO NOT MODIFY - This file is automatically generated during compilation.

class WaveParameters : public Parameters {
  public:
    enum InitialConditions{
      GAUSSIAN,
      FLAT,
    };

    WaveParameters() : Parameters(1){
      mInitialConditions = GAUSSIAN;
      mKOSigma = 0.0;
      mGaussianAmplitude = 1.0;
    }

    inline void setInitialConditions(InitialConditions val){
      mInitialConditions = val;
    }

    inline InitialConditions getInitialConditions(){
      return mInitialConditions;
    }

    inline void setKOSigma(double KOSigma){
      mKOSigma = KOSigma;
    }

    inline double getKOSigma(){
      return mKOSigma;
    }

    inline void setGaussianAmplitude(double GaussianAmplitude){
      mGaussianAmplitude = GaussianAmplitude;
    }

    inline double getGaussianAmplitude(){
      return mGaussianAmplitude;
    }

  private:
    InitialConditions mInitialConditions;
    double mKOSigma;
    double mGaussianAmplitude;
};

#endif

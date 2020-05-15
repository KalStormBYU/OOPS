#include <waveparser.h>
#include <iostream>

using namespace std;

void WaveParser::updateParameters(string fname, Parameters* params){
    if(!checkId(params)){
        return;
    }
    reader.clearData();
    ParamReader::ParamResult result = reader.readFile(fname);
    if(result != ParamReader::SUCCESS){
        cout << "An error occurred while trying to read " << fname << ".\n";
        return;
    }
    WaveParameters *pars = (WaveParameters*) params;

    if(!reader.hasSection(string("Wave"))){
        return;
    }

    if(reader.hasParameter(string("Wave"),string("npoints"))){
        int result = reader.readAsInt(string("Wave"),string("npoints"));
        if(result <= 50001 && result >= 5){
            pars->setnpoints(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("xmin"))){
        double result = reader.readAsDouble(string("Wave"),string("xmin"));
        if(result <= 1.000000e+03 && result >= -1.000000e+03){
            pars->setxmin(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("xmax"))){
        double result = reader.readAsDouble(string("Wave"),string("xmax"));
        if(result <= 1.000000e+03 && result >= -1.000000e+03){
            pars->setxmax(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("tmax"))){
        double result = reader.readAsDouble(string("Wave"),string("tmax"));
        if(result <= 1.000000e+04 && result >= 0.000000e+00){
            pars->settmax(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("cfl"))){
        double result = reader.readAsDouble(string("Wave"),string("cfl"));
        if(result <= 1.000000e+00 && result >= 0.000000e+00){
            pars->setcfl(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("output_frequency"))){
        int result = reader.readAsInt(string("Wave"),string("output_frequency"));
        if(result <= 1000000 && result >= 1){
            pars->setoutput_frequency(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("id_gauss_amp"))){
        double result = reader.readAsDouble(string("Wave"),string("id_gauss_amp"));
        if(result <= 1.000000e+01 && result >= 0.000000e+00){
            pars->setid_gauss_amp(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("id_gauss_center"))){
        double result = reader.readAsDouble(string("Wave"),string("id_gauss_center"));
        if(result <= 1.000000e+03 && result >= -1.000000e+03){
            pars->setid_gauss_center(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("id_gauss_width"))){
        double result = reader.readAsDouble(string("Wave"),string("id_gauss_width"));
        if(result <= 1.000000e+01 && result >= 1.000000e-03){
            pars->setid_gauss_width(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

    if(reader.hasParameter(string("Wave"),string("KOSigma"))){
        double result = reader.readAsDouble(string("Wave"),string("KOSigma"));
        if(result <= 1.000000e+00 && result >= 0.000000e+00){
            pars->setKOSigma(result);
        }
        else{
            cout << "Parameter %s out of range.\n";
        }
    }

}

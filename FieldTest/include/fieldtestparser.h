#ifndef FIELDTEST_PARSER_H
#define FIELDTEST_PARSER_H

#include <paramparser.h>
#include <fieldtestparameters.h>

// DO NOT MODIFY - This file is automatically generated during compilation.

class FieldTestParser : ParamParser{
  public:
    FieldTestParser() : ParamParser(1){}

    virtual ~FieldTestParser(){}

    virtual void updateParameters(std::string fname, Parameters* params);
};

#endif

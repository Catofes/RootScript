//
// Created by herbertqiao on 11/17/16.
//

#ifndef ROOTSCRIPT_FUNCTION_H
#define ROOTSCRIPT_FUNCTION_H

#include <iostream>
#include <string.h>
#include <vector>
#include <RooRealVar.h>
#include <RooAbsPdf.h>

using namespace std;

class Function
{
public:
    Function()
    {};

    ~Function();

    Function(const string &name);

    void SetName(const string &name)
    { _name = name; }

    string GetName()
    { return _name; }

    void AddParameter(const string &name, const double &mean, const double &width, const bool &adaptable = true);

    bool AdaptParameter();

    virtual bool Create()
    { return true; }

    virtual string Print()
    { return ""; }

protected:
    string _name = "";
    vector<RooRealVar> _parameter_list;
    vector<bool> _parameter_list_adaptable;
    RooAbsPdf *_pdf = NULL;
};


#endif //ROOTSCRIPT_FUNCTION_H

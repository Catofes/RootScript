//
// Created by herbertqiao on 11/17/16.
//

#include "Function.h"

Function::Function(const string &name)
        : _name(name)
{}

Function::~Function()
{
    if (_pdf != NULL)
        delete _pdf;
}

void Function::AddParameter(const string &name, const double &mean, const double &width, const bool &adaptable)
{
    string full_name = _name + "_" + name;
    RooRealVar var(full_name.c_str(), full_name.c_str(), mean, mean - width / 2, mean + width / 2);
    bool adaptable_ = adaptable;
    _parameter_list.push_back(var);
    _parameter_list_adaptable.push_back(adaptable_);
}

bool Function::AdaptParameter()
{
    bool adapted = false;
    for (int i = 0; i < _parameter_list.size(); i++) {
        if _parameter_list_adaptable
    }

}
//
// Created by herbertqiao on 12/13/16.
//

#include "G4DataConvert.h"

G4DataConvert::G4DataConvert(const string &path)
{
    chain = new TChain("mcTree");
    chain->Add(path.c_str());
    chain->SetBranchAddress("xd", &x);
    chain->SetBranchAddress("yd", &y);
    chain->SetBranchAddress("zd", &z);
    chain->SetBranchAddress("energy", &e);
}

G4DataConvert::~G4DataConvert()
{
    if (chain != NULL) {
        delete chain;
    }
}

int G4DataConvert::get_total_entries()
{
    return chain->GetEntries();
}

vector<hit> G4DataConvert::get_hits(int entry)
{
    vector<hit> result;
    if (entry >= 0 && entry < get_total_entries())
        chain->GetEntry(entry);
    else
        throw std::out_of_range("Entry out of range.");
    for (int i = 0; i < x->size(); i++)
        result.push_back(make_tuple((*x)[i], (*y)[i], (*z)[i], (*e)[i]));
    return result;
}

vector<double> *G4DataConvert::get_hit_x(int entry)
{
    if (entry >= 0 && entry < get_total_entries())
        chain->GetEntry(entry);
    else
        throw std::out_of_range("Entry out of range.");
    return x;
}

vector<double> *G4DataConvert::get_hit_y(int entry)
{
    if (entry >= 0 && entry < get_total_entries())
        chain->GetEntry(entry);
    else
        throw std::out_of_range("Entry out of range.");
    return y;
}

vector<double> *G4DataConvert::get_hit_z(int entry)
{
    if (entry >= 0 && entry < get_total_entries())
        chain->GetEntry(entry);
    else
        throw std::out_of_range("Entry out of range.");
    return z;
}

vector<double> *G4DataConvert::get_hit_e(int entry)
{
    if (entry >= 0 && entry < get_total_entries())
        chain->GetEntry(entry);
    else
        throw std::out_of_range("Entry out of range.");
    return e;
}

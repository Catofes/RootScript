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
    if(chain != NULL){
        
    }
}
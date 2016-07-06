//
// Created by herbertqiao on 6/28/16.
//

#include "Raw2Info.hh"
#include "TTree.h"

int main(char *args, char **argv)
{
    Raw2Info r(argv[0]);

}

Raw2Info::Raw2Info(string input)
{
    file = new TFile(input.c_str(), "Read");
    if (!file)
        throw runtime_error("Error File");
    tree = (TTree *) file->Get("mcTree");
    if (!tree)
        throw runtime_error("Error Tree");
    entries = tree->GetEntries();
    tree->SetBranchAddress("energy", &_energy);
    tree->SetBranchAddress("xd", &_xd);
    tree->SetBranchAddress("yd", &_yd);
    tree->SetBranchAddress("zd", &_zd);
    tree->SetBranchAddress("td", &_td);
    tree->SetBranchAddress("volume", &_volume);
    tree->SetBranchAddress("primaryType", &_primaryType);
}

void Raw2Info::Process(int i)
{
    if (i >= entries)
        throw runtime_error("Out ot border");
    tree->GetEntry(i);

}
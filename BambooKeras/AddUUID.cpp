//
// Created by herbertqiao on 3/13/17.
//

#include <TFile.h>
#include <TTree.h>
#include "sole.hpp"
#include "iostream"

using namespace std;

int main(int argc, char *argv[])
{
    if (argc < 2)
        throw std::runtime_error("Missing parameters.");
    TFile *file = new TFile(argv[1], "Update");
    TTree *tree = (TTree *) file->Get("mcTree");
    string uuid = "";
    tree->Branch("uuid", &uuid);
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        uuid = sole::uuid4().str();
        tree->GetBranch("uuid")->Fill();
    }
    tree->Write();
    file->Close();
}
//
// Created by herbertqiao on 11/30/16.
//

#include "LengthDistribution.h"
#include <iostream>
#include <TChain.h>
#include <TH1F.h>
#include <TFile.h>

using namespace std;


void LengthDistribution(const char *input_, const char *output_)
{
    TChain *chain = new TChain("mcTree");
    chain->Add(input_);
    long entries = chain->GetEntries();
    cout << "Total Entries:" << entries << endl;
    vector<double> *xd = 0, *yd = 0, *zd = 0;
    double x, y, z, energy;
    chain->SetBranchAddress("xd", &xd);
    chain->SetBranchAddress("yd", &yd);
    chain->SetBranchAddress("zd", &zd);
    chain->SetBranchAddress("totalEnergy", &energy);

    TFile *f = new TFile(output_, "RECREATE");
    TTree *t = new TTree("length_distribution", "Length Distribution");
    t->Branch("x", &x);
    t->Branch("y", &y);
    t->Branch("z", &z);
    t->Branch("energy", &energy);

    for (int i = 0; i < entries; i++) {
        if (i % 10000 == 0)
            cout << i << "/" << entries << endl;
        chain->GetEntry(i);
        long n_hits = xd->size();
        if (n_hits <= 0)
            continue;
        double x_min, x_max, y_min, y_max, z_min, z_max;
        x_min = (*xd)[0];
        x_max = (*xd)[0];
        y_min = (*yd)[0];
        y_max = (*yd)[0];
        z_min = (*zd)[0];
        z_max = (*zd)[0];
        for (int j = 0; j < n_hits; j++) {
            if ((*xd)[j] > x_max)
                x_max = (*xd)[j];
            if ((*xd)[j] < x_min)
                x_min = (*xd)[j];
            if ((*yd)[j] > y_max)
                y_max = (*yd)[j];
            if ((*yd)[j] < y_min)
                y_min = (*yd)[j];
            if ((*zd)[j] > z_max)
                z_max = (*zd)[j];
            if ((*zd)[j] < z_min)
                z_min = (*zd)[j];
        }
        x = x_max - x_min;
        y = y_max - y_min;
        z = z_max - z_min;
        t->Fill();
    }
    t->Write();
    f->Close();
}


int main(int argc, char *argv[])
{
    if (argc != 3)
        throw std::runtime_error("Error Input. Use LengthDistribution <input file> <output file>");
    cout << "Input File:" << argv[1] << endl;
    cout << "Output File:" << argv[2] << endl;
    LengthDistribution(argv[1], argv[2]);
}
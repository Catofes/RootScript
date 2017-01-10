#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <iostream>

using namespace std;

int main()
{
    TFile *f = new TFile("input_data.root", "RECREATE");
    TTree *t = new TTree("Data", "Data");
    vector<double> x, y, z, weight;
    cout << "a";
    t->Branch("x", &x);
    t->Branch("y", &y);
    t->Branch("z", &z);
    t->Branch("weight", &weight);
    TRandom random;
    for (int j = 0; j < 100; j++) {
        x.clear();
        y.clear();
        z.clear();
        weight.clear();
        for (int i = 0; i < 100; i++) {
            x.push_back(random.Rndm());
            y.push_back(random.Rndm());
            z.push_back(random.Rndm());
            weight.push_back(random.Rndm());
        }
        t->Fill();
    }
    f->Write();
}
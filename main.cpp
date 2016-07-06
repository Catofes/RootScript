#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TTree.h>
#include <TCanvas.h>

using namespace std;

int main()
{
    TChain *t = new TChain("simple_out");
    t->Add("simple_output.root");
    TCanvas *canvas = new TCanvas("Energy in keV", "Energy in keV", 600, 400);
    TH1F *h1 = new TH1F("Energy in keV", "Energy in keV", 100, 2457 - 200, 2457 + 200);
    double totally_energy;
    t->SetBranchAddress("energy", &totally_energy);
    int nEntries = t->GetEntries();
    cout << nEntries << endl;
    for (int i = 0; i < nEntries; i++) {
        t->GetEntry(i);
        if (totally_energy > 2457 - 200 && totally_energy < 2457 + 200) {
            h1->Fill(totally_energy);
        }
    }
    h1->Draw();
    canvas->Update();
    TLine *l = new TLine(2457.83 - 41.75, 0, 2457.83 - 41.75, canvas->GetUymax());
    l->SetLineColor(kRed);
    l->Draw();
    TLine *l = new TLine(2457.83 + 41.75, 0, 2457.83 + 41.75, canvas->GetUymax());
    l->SetLineColor(kRed);
    l->Draw();
}
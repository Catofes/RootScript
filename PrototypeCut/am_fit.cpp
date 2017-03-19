#include <iostream>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>

using namespace std;

void fit(const string &input_file)
{
    TChain *chain = new TChain("SmearEnergy");
    chain->Add(input_file.c_str());
    TH1F *h = new TH1F("energy", "Energy", 500, 0, 70);
    chain->Draw("smearEnergy>>energy", "smearEnergy>0&&primaryType==\"Am241\"");
    double par[9];
    TF1 *g1 = new TF1("g1", "gaus", 12, 16);
    TF1 *g2 = new TF1("g2", "gaus", 16, 20);
    TF1 *total = new TF1("total", "gaus(0)+gaus(3)", 12, 20);
    total->SetLineColor(2);
    h->Fit(g1, "R0");
    h->Fit(g2, "R+");
    g1->GetParameters(&par[0]);
    g2->GetParameters(&par[3]);
    total->SetParameters(par);
    h->Fit(total, "R");
    cout << "Integral(12,20): " << total->Integral(12, 20) << endl;
    TF1 *g3 = new TF1("g1", "gaus", 50, 70);
    g3->SetLineColor(3);
    h->Fit(g3, "R+");
    cout << "Integral(50,70): " << g3->Integral(50, 70) << endl;
    TF1 *g4 = new TF1("g2", "gaus", 25, 35);
    g4->SetLineColor(4);
    h->Fit(g4, "R+");
    cout << "Integral(50,70): " << g4->Integral(25, 35) << endl;
}

int main(int argc, char *argv[])
{
    fit(string(argv[1]));
}
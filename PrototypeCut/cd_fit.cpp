#include <iostream>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>
#include <TApplication.h>
#include <argparse.h>

using namespace std;

void fit(const string &input_file)
{
    TChain *chain = new TChain("SmearEnergy");
    chain->Add(input_file.c_str());
    TH1F *h = new TH1F("energy", "Energy", 500, 0, 70);
    chain->Draw("smearEnergy>>energy", "smearEnergy>0");
    TF1 *g1 = new TF1("g1", "gaus", 4, 12);
    TF1 *g2 = new TF1("g2", "gaus", 14, 32);
    TF1 *g3 = new TF1("g3", "gaus", 36, 54);
    g1->SetLineColor(2);
    g2->SetLineColor(3);
    g3->SetLineColor(4);
    h->Fit(g1, "R");
    cout << "Integral(4,12): " << g1->Integral(4, 12) << endl;
    h->Fit(g2, "R+");
    cout << "Integral(14,32): " << g2->Integral(14, 32) << endl;
    h->Fit(g3, "R+");
    cout << "Integral(36,54): " << g3->Integral(36, 54) << endl;
}

int main(int argc, char *argv[])
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.parse(argc, argv);
    TApplication *myapp = new TApplication("App", &argc, argv);
    fit(parser.retrieve<string>("i"));
    myapp->Run();
}
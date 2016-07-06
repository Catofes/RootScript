//
// Created by herbertqiao on 7/6/16.
//

#include <TH1.h>
#include <TRandom3.h>
#include "SmearMC.hh"

using namespace std;

SmearMC::SmearMC(std::string input_path, std::string output_path)
{
    _input_chain = new TChain("t1");
    _input_chain->Add(input_path.c_str());
    _output_file = new TFile(output_path.c_str(), "RECREATE");
    _output_tree = new TTree("t1", "Energy results");
    _input_chain->SetBranchAddress("Cha_Energy", &_cha_energy);
    _output_tree->Branch("Cha_Energy", &_cha_energy, "Cha_Energy/F");
}

void SmearMC::smear()
{
    TH1F *energy = new TH1F("he", "he", 6000, 0.0, 3000.0);
    energy->SetXTitle("E [keV]");
    energy->SetYTitle("Entries / 0.5keV");
    TH1F *hparameter = new TH1F("hparameter", "hparameter", 12, 0.0, 12.0);
    int entries_num = _input_chain->GetEntries();
    cout << "Events Number:" << entries_num;

    Float_t sigma_e;
    TRandom3 r;
    for (int i = 0; i < entries_num; i++) {
        _input_chain->GetEntry(i);
        sigma_e = 0.512 + 0.0002632 * _cha_energy;
        _cha_energy += r.Gaus(0.0, sigma_e);
        energy->Fill(_cha_energy);
        _output_tree->Fill();
    }
    hparameter->SetBinContent(1, _simulate_num / 1000000);//runtime
    hparameter->SetBinContent(2, _simulate_num / 1000000);//livetime
    _output_file->cd();
    energy->Write();
    hparameter->Write();
    _output_tree->Write();
    _output_file->Close();
}

int main(int argc, char *argv[])
{
    if (argc < 2)
        throw std::runtime_error("Error Input. Use MergeMC input_file output_file.");
    SmearMC(argv[0], argv[1]).smear();
}
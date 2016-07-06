//
// Created by herbertqiao on 7/6/16.
//

#include <TH1.h>
#include "MergeMC.hh"

using namespace std;

MergeMC::MergeMC(std::string input_path, std::string output_path)
{
    _input_chain = new TChain("mcTree");
    _input_chain->Add(input_path.c_str());
    _output_file = new TFile(output_path.c_str(), "RECREATE");
    _output_tree = new TTree("t1", "Energy results");
    _input_chain->SetBranchAddress("totalEnergy", &_total_energy);
    _output_tree->Branch("Cha_Energy", &_cha_energy, "Cha_Energy");
}

void MergeMC::merge()
{
    _output_file->cd();
    TH1F *energy = new TH1F("he", "he", 6000, 0.0, 3000.0);
    energy->SetXTitle("E [keV]");
    energy->SetYTitle("Entries / 0.5keV");
    TH1F *hparameter = new TH1F("hparameter", "hparameter", 12, 0.0, 12.0);
    int entries_num = _input_chain->GetEntries();
    cout << "Events Number:" << entries_num;
    if (entries_num == 0)
        throw std::runtime_error("Zero Events");
    for (int i = 0; i < entries_num; i++) {
        _input_chain->GetEntry(i);
        _cha_energy = _total_energy;
        energy->Fill(_total_energy);
        _output_tree->Fill();
    }
    hparameter->SetBinContent(1, _simulate_num / 1000000);//runtime
    hparameter->SetBinContent(2, _simulate_num / 1000000);//livetime
    energy->Write();
    hparameter->Write();
    _output_tree->Write();
    _output_file->Close();
}

int main(int argc, char *argv[])
{
    if (argc != 3)
        throw std::runtime_error("Error Input. Use MergeMC input_file output_file.");
    cout << "Input File:" << argv[1] << endl;
    cout << "Output File:" << argv[2] << endl;
    MergeMC(argv[1], argv[2]).merge();
}
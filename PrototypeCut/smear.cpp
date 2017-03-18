//
// Created by herbertqiao on 2/13/17.
//

#include "base_convert.h"
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>

class SmearTotalEnergy
        : public BaseConvert
{
public:
    SmearTotalEnergy(const string &input_path, const double &resolution, const string &output_path);

    void process(int i);

    void final();

    double energy_in_readout();

private:
    TFile *_out_file = 0;
    TTree *_tree = 0;
    double _smear_energy = 0;
    double _readout_energy = 0;
    double _resolution = 0;
    TRandom *_rand = new TRandom();
};

SmearTotalEnergy::SmearTotalEnergy(const string &input_path, const double &resolution = 0.03,
                                   const string &output_path = "output.root")
        : BaseConvert(input_path)
{
    _out_file = new TFile(output_path.c_str(), "RECREATE");
    _tree = new TTree("SmearEnergy", "SmearEnergy");
    _tree->Branch("totalEnergy", &totalEnergy);
    _tree->Branch("readoutEnergy", &_readout_energy);
    _tree->Branch("smearEnergy", &_smear_energy);
    _tree->Branch("primaryType", &primaryType);
    _resolution = resolution;
}

double SmearTotalEnergy::energy_in_readout()
{
    double total_energy = 0;
    for (int i = 0; i < xd->size(); i++) {
        double x = (*xd)[i];
        double y = (*zd)[i];
        double z = (*yd)[i];
        double e = (*energy)[i];
        if (-100 < x && x < 100 && -100 < y && y < 100)
            total_energy += e;
    }
    return total_energy;
}

void SmearTotalEnergy::process(int entry)
{
    chain->GetEntry(entry);
    _readout_energy = energy_in_readout();
    if (_readout_energy > 0)
        _smear_energy = _rand->Gaus(_readout_energy,
                                    _readout_energy * _resolution / 2.355 * sqrt(2457.83 / _readout_energy));
    if (_smear_energy < 0)
        _smear_energy = 0;
    _tree->Fill();
}

void SmearTotalEnergy::final()
{
    process_all();
    _out_file->cd();
    _tree->Write();
    _out_file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-r", "--resolution", 1, false);
    parser.addArgument("-o", "--output", 1, false);
    parser.parse(argc, argv);
    //TApplication *myapp = new TApplication("App", &argc, argv);
    SmearTotalEnergy convert(parser.retrieve<string>("input"), stod(parser.retrieve<string>("resolution")),
                             parser.retrieve<string>("output"));
    cout << "Bind finished" << endl;
    convert.final();
    //myapp->Run();
}
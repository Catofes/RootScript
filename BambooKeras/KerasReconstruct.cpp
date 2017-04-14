//
// Created by herbertqiao on 3/13/17.
//

#include <TChain.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <argparse.h>
#include "ReadoutWave.hh"
#include "ReadoutStruct.hh"
#include "base_convert.h"

using namespace std;

class ReconstructConvert
        : public BaseConvert
{
public:
    ReconstructConvert(const string &input_file_name, const string &output_file_name,
                       const string &json_file_name = "main.json");

    void final();

    void process(int i);

private:
    TFile *output_file;
    TTree *output_tree;
    ReadoutWave *_readout_wave = 0;
    double _trigger_Energy;
    bool _triggered;
    bool _cross_cathode;
    int _gap_count;
    int _out_pixel_count;
    int _record_count;

    double work_function = 21.9 / 1000;
    double drift_velocity = 1.87139;

    vector<double> *_xzx = new vector<double>;
    vector<double> *_xzz = new vector<double>;
    vector<double> *_xze = new vector<double>;
    vector<double> *_yzy = new vector<double>;
    vector<double> *_yzz = new vector<double>;
    vector<double> *_yze = new vector<double>;
    ReadoutPlane *plane = 0;

};

ReconstructConvert::ReconstructConvert(const string &input_file_name, const string &output_file_name,
                                       const string &json_file_name)
        : BaseConvert(input_file_name, "ReadoutWave")
{
    chain->SetBranchAddress("triggerEnergy", &_trigger_Energy);
    chain->SetBranchAddress("readoutWave", &_readout_wave);
    chain->SetBranchAddress("triggered", &_triggered);
    chain->SetBranchAddress("gapCount", &_gap_count);
    chain->SetBranchAddress("outPixelCount", &_out_pixel_count);
    chain->SetBranchAddress("recordCount", &_record_count);
    chain->SetBranchAddress("crossCathode", &_cross_cathode);

    output_file = new TFile(output_file_name.c_str(), "RECREATE");
    output_tree = new TTree("MLData", "MachineLearning Data");

    output_tree->Branch("runId", &runId, "runId/I");
    output_tree->Branch("eventId", &eventId, "eventId/I");
    output_tree->Branch("nHits", &nHits, "nHits/I");
    output_tree->Branch("trackId", &trackId);
    output_tree->Branch("parentId", &parentId);
    output_tree->Branch("type", &type);
    output_tree->Branch("parentType", &parentType);
    output_tree->Branch("creatorProcess", &creatorProcess);
    output_tree->Branch("depositionProcess", &depositionProcess);
    output_tree->Branch("volume", &volume);
    output_tree->Branch("totalEnergy", &totalEnergy);
    output_tree->Branch("xd", &xd);
    output_tree->Branch("yd", &yd);
    output_tree->Branch("zd", &zd);
    output_tree->Branch("td", &td);
    output_tree->Branch("energy", &energy);
    output_tree->Branch("nPrimaries", &nPrimaries);
    output_tree->Branch("primaryType", &primaryType);
    output_tree->Branch("primaryId", &primaryId);
    output_tree->Branch("primaryEnergy", &primaryEnergy);
    output_tree->Branch("primaryPx", &primaryPx);
    output_tree->Branch("primaryPy", &primaryPy);
    output_tree->Branch("primaryPz", &primaryPz);
    output_tree->Branch("primaryX", &primaryX);
    output_tree->Branch("primaryY", &primaryY);
    output_tree->Branch("primaryZ", &primaryZ);
    output_tree->Branch("uuid", &uuid);

    output_tree->Branch("triggerEnergy", &_trigger_Energy);
    output_tree->Branch("readoutWave", &_readout_wave);
    output_tree->Branch("triggered", &_triggered);
    output_tree->Branch("gapCount", &_gap_count);
    output_tree->Branch("outPixelCount", &_out_pixel_count);
    output_tree->Branch("recordCount", &_record_count);
    output_tree->Branch("crossCathode", &_cross_cathode);

    output_tree->Branch("xzx", &_xzx);
    output_tree->Branch("xzz", &_xzz);
    output_tree->Branch("xze", &_xze);
    output_tree->Branch("yzy", &_yzy);
    output_tree->Branch("yzz", &_yzz);
    output_tree->Branch("yze", &_yze);
    plane = new ReadoutPlane(json_file_name);
}

void ReconstructConvert::process(int i)
{
    chain->GetEntry(i);
    if (_trigger_Energy < 2400)
        return;
    _xzx->clear();
    _xzz->clear();
    _xze->clear();
    _yzy->clear();
    _yzz->clear();
    _yze->clear();
    for (const auto &u:(*_readout_wave).detectors) {
        int channel_id = u.first;
        if ((channel_id % 1000) < 100) {
            double x = plane->GetX(channel_id);
            for (int j = 0; j < u.second.size(); j++) {
                if (u.second[j] <= 0)
                    continue;
                double e = u.second[j] * work_function;
                double z = j * 0.2 * drift_velocity;
                _xzx->push_back(x);
                _xzz->push_back(z);
                _xze->push_back(e);
            }
        } else {
            double y = plane->GetY(channel_id);
            for (int j = 0; j < u.second.size(); j++) {
                if (u.second[j] <= 0)
                    continue;
                double e = u.second[j] * work_function;
                double z = j * 0.2 * drift_velocity;
                _yzy->push_back(y);
                _yzz->push_back(z);
                _yze->push_back(e);
            }
        }
    }
    output_tree->Fill();
}

void ReconstructConvert::final()
{
    process_all();
    output_tree->Write();
    output_file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-j", "--json", 1, false);
    parser.addArgument("-o", "--output", 1, false);
    parser.parse(argc, argv);
    ReconstructConvert(parser.retrieve<string>("input"), parser.retrieve<string>("output"),
                       parser.retrieve<string>("json")).final();
}
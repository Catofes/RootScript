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

using namespace std;

void convert(const string &input_file_name, const string &output_file_name, const string &json_file_name = "main.json")
{
    double work_function = 21.9 / 1000;
    double drift_velocity = 1.87139;

    TChain *chain = new TChain("ReadoutWave");
    chain->Add(input_file_name.c_str());
    ReadoutWave *_readout_wave = 0;
    int _runId = 0;
    int _eventId = 0;
    double _totalEnergy = 0;
    string *_uuid = 0;
    double _trigger_Energy = 0;
    bool _triggered = 0;
    bool _cross_cathode = 0;
    int _gap_count = 0;
    int _out_pixel_count = 0;
    int _record_count = 0;

    chain->SetBranchAddress("runId", &_runId);
    chain->SetBranchAddress("eventId", &_eventId);
    chain->SetBranchAddress("totalEnergy", &_totalEnergy);
    chain->SetBranchAddress("uuid", &_uuid);
    chain->SetBranchAddress("triggerEnergy", &_trigger_Energy);
    chain->SetBranchAddress("readoutWave", &_readout_wave);
    chain->SetBranchAddress("triggered", &_triggered);
    chain->SetBranchAddress("gapCount", &_gap_count);
    chain->SetBranchAddress("outPixelCount", &_out_pixel_count);
    chain->SetBranchAddress("recordCount", &_record_count);
    chain->SetBranchAddress("crossCathode", &_cross_cathode);

    vector<double> *_xzx = new vector<double>;
    vector<double> *_xzz = new vector<double>;
    vector<double> *_xze = new vector<double>;
    vector<double> *_yzy = new vector<double>;
    vector<double> *_yzz = new vector<double>;
    vector<double> *_yze = new vector<double>;

    TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
    TTree *output_tree = new TTree("MLData", "MachineLearning Data");
    output_tree->Branch("xzx", &_xzx);
    output_tree->Branch("xzz", &_xzz);
    output_tree->Branch("xze", &_xze);
    output_tree->Branch("yzy", &_yzy);
    output_tree->Branch("yzz", &_yzz);
    output_tree->Branch("yze", &_yze);
    output_tree->Branch("runId", &_runId, "runId/I");
    output_tree->Branch("eventId", &_eventId, "eventId/I");
    output_tree->Branch("totalEnergy", &_totalEnergy, "totalEnergy/D");
    output_tree->Branch("uuid", &_uuid);
    output_tree->Branch("triggerEnergy", &_trigger_Energy, "triggerEnergy/D");
    output_tree->Branch("triggered", &_triggered);
    output_tree->Branch("gapCount", &_gap_count);
    output_tree->Branch("outPixelCount", &_out_pixel_count);
    output_tree->Branch("recordCount", &_record_count);
    output_tree->Branch("crossCathode", &_cross_cathode);


    ReadoutPlane *plane = new ReadoutPlane(json_file_name);
    for (int i = 0; i < chain->GetEntries(); i++) {
        if (i % 10 == 0) {
            cout << i << "/" << chain->GetEntries() << endl;
        }
        chain->GetEntry(i);
        if (_trigger_Energy < 2400)
            continue;
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
    convert(parser.retrieve<string>("input"), parser.retrieve<string>("output"), parser.retrieve<string>("json"));
}
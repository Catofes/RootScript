//
// Created by herbertqiao on 3/16/17.
//

#include <iostream>
#include <TChain.h>
#include <vector>
#include <map>
#include <string>
#include <TFile.h>
#include <fstream>
#include <argparse.h>

using namespace std;

const double cluster_x_size = 3;
const double cluster_z_size = 2.99;
const double energy_cut_low = 2436.9;
const double energy_cut_high = 2478.7;

pair<double, double> get_center(const vector<double> &x, const vector<double> &z, const vector<double> &e)
{
    double total_energy = 0;
    double total_x = 0;
    double total_z = 0;
    for (int i = 0; i < x.size(); i++) {
        total_energy += e[i];
        total_x += x[i] * e[i];
        total_z += z[i] * e[i];
    }
    return make_pair(total_x / total_energy, total_z / total_energy);
};

double cluster(map<pair<int, int>, double> &cluster_data, const vector<double> &x, const vector<double> &z,
               const vector<double> &e)
{
    double center_x = 0, center_z = 0;
    tie(center_x, center_z) = get_center(x, z, e);
    for (int i = 0; i < x.size(); i++) {
        auto cluster_x = (int) floor((x[i] - center_x) / cluster_x_size);
        auto cluster_z = (int) floor((z[i] - center_z) / cluster_z_size);
        auto cluster_id = make_pair(cluster_x, cluster_z);
        if (cluster_data.find(cluster_id) != cluster_data.end())
            cluster_data[cluster_id] += e[i];
        else
            cluster_data.insert(make_pair(cluster_id, e[i]));
    }
    double max_energy = 0;
    for (auto &u:cluster_data) {
        if (u.second > max_energy)
            max_energy = u.second;
    }

    return max_energy;
}

void convert(const string &input_file_name, const string &output_file_name, const string &type)
{
    TChain *chain = new TChain("MLData");
    chain->Add(input_file_name.c_str());

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
    vector<double> *_xzx = 0;
    vector<double> *_xzz = 0;
    vector<double> *_xze = 0;
    vector<double> *_yzy = 0;
    vector<double> *_yzz = 0;
    vector<double> *_yze = 0;

    chain->SetBranchAddress("xzx", &_xzx);
    chain->SetBranchAddress("xzz", &_xzz);
    chain->SetBranchAddress("xze", &_xze);
    chain->SetBranchAddress("yzy", &_yzy);
    chain->SetBranchAddress("yzz", &_yzz);
    chain->SetBranchAddress("yze", &_yze);
    chain->SetBranchAddress("runId", &_runId);
    chain->SetBranchAddress("eventId", &_eventId);
    chain->SetBranchAddress("totalEnergy", &_totalEnergy);
    chain->SetBranchAddress("uuid", &_uuid);
    chain->SetBranchAddress("triggerEnergy", &_trigger_Energy);
    chain->SetBranchAddress("triggered", &_triggered);
    chain->SetBranchAddress("gapCount", &_gap_count);
    chain->SetBranchAddress("outPixelCount", &_out_pixel_count);
    chain->SetBranchAddress("recordCount", &_record_count);
    chain->SetBranchAddress("crossCathode", &_cross_cathode);

    vector<int> *_cluster_xzx = new vector<int>;
    vector<int> *_cluster_xzz = new vector<int>;
    vector<double> *_cluster_xze = new vector<double>;
    vector<int> *_cluster_yzy = new vector<int>;
    vector<int> *_cluster_yzz = new vector<int>;
    vector<double> *_cluster_yze = new vector<double>;
    string _type = type;
    double _cluster_max_energy = 0;
    TFile *output_file = new TFile(output_file_name.c_str(), "RECREATE");
    TTree *output_tree = new TTree("PictureData", "Picture Data");
    output_tree->Branch("cluster_xzx", &_cluster_xzx);
    output_tree->Branch("cluster_xzz", &_cluster_xzz);
    output_tree->Branch("cluster_xze", &_cluster_xze);
    output_tree->Branch("cluster_yzy", &_cluster_yzy);
    output_tree->Branch("cluster_yzz", &_cluster_yzz);
    output_tree->Branch("cluster_yze", &_cluster_yze);
    output_tree->Branch("cluster_max_energy", &_cluster_max_energy);
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
    output_tree->Branch("type", &_type);


    for (int i = 0; i < chain->GetEntries(); i++) {
        if (i % 1000 == 0)
            cout << i << "/" << chain->GetEntries() << endl;
        chain->GetEntry(i);
        if (_trigger_Energy < energy_cut_low || _trigger_Energy > energy_cut_high)
            continue;
        map<pair<int, int>, double> cluster_info;
        _cluster_xzx->clear();
        _cluster_xzz->clear();
        _cluster_xze->clear();
        _cluster_yzy->clear();
        _cluster_yzz->clear();
        _cluster_yze->clear();
        auto cluster_xz_energy = cluster(cluster_info, (*_xzx), (*_xzz), (*_xze));
        for (auto &u:cluster_info) {
            _cluster_xzx->push_back(u.first.first);
            _cluster_xzz->push_back(u.first.second);
            _cluster_xze->push_back(u.second);
        }
        cluster_info.clear();
        auto cluster_yz_energy = cluster(cluster_info, (*_yzy), (*_yzz), (*_yze));
        for (auto &u:cluster_info) {
            _cluster_yzy->push_back(u.first.first);
            _cluster_yzz->push_back(u.first.second);
            _cluster_yze->push_back(u.second);
        }
        _cluster_max_energy = max(cluster_xz_energy, cluster_yz_energy);
        output_tree->Fill();
    }
    output_tree->Write();
    output_file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-o", "--output", 1, false);
    parser.addArgument("-t", "--type", 1, false);
    parser.parse(argc, argv);
    convert(parser.retrieve<string>("input"), parser.retrieve<string>("output"), parser.retrieve<string>("type"));
}
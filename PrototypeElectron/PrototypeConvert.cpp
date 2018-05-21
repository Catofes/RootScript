//
// Created by herbertqiao on 4/19/17.
//

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <TRandom.h>
#include <TMath.h>
#include "base_convert.h"
#include "json.h"
#include "TFile.h"

using namespace std;

typedef tuple<double, double> MicroMegas; //x,y in mm

class PrototypeConvert
        : public BaseConvert
{
public:
    PrototypeConvert(string &input_path, string &output_path, const string &json_path);

    void load_parameter(string &input_path);

    void final();

    void process(int i);

private:
    int get_electron_numbers(double energy);

    bool out_of_border(double x, double y);

    void get_electron_info();

    tuple<int, int, int> get_trigger_info(vector<double> &electron_info);

    TFile *_output_file;
    TTree *_output_tree;
    TRandom random;

    double _work_function;
    double _fano;
    double _transverse_diffusion;
    double _longitudinal_diffusion;
    double _drift_velocity;
    double _z_plane;
    double _trigger_threshold;
    double _micromegas_size;
    double _trigger_bin_length;
    double _trigger_bin_size;
    double _trigger_bin_offset;

    int _gap_count;
    int _record_count;
    bool _triggered;
    double _trigger_energy;

    vector<double> _electron_info;
    vector<MicroMegas> _micromegas_info;
};

PrototypeConvert::PrototypeConvert(string &input_path, string &json_path, const string &output_path)
        : BaseConvert(input_path, "mcTree"), random(), _electron_info()
{
    _work_function = 21.9 / 1000;
    _fano = 0.14;
    _transverse_diffusion = 0.0101968; // cm1/2
    _longitudinal_diffusion = 0.0139049; // cm1/2
    _drift_velocity = 1.87139; // mm/us
    _z_plane = 392.5; //mm
    _trigger_threshold = 10; //keV
    _trigger_bin_length = 0.2; //ms
    _trigger_bin_size = 512;
    _trigger_bin_offset = 256;
    _micromegas_size = 200; //mm

    load_parameter(json_path);

    _output_file = new TFile(output_path.c_str(), "RECREATE");
    _output_tree = new TTree("DriftResult", "DriftResult");

    _output_tree->Branch("totalEnergy", &totalEnergy);
    _output_tree->Branch("triggerEnergy", &_trigger_energy, "triggerEnergy/D");
    _output_tree->Branch("triggered", &_triggered);
    _output_tree->Branch("gapCount", &_gap_count);
    _output_tree->Branch("recordCount", &_record_count);
    _output_tree->Branch("primaryType", &primaryType);
}

void PrototypeConvert::load_parameter(string &input_path)
{
    Json::Reader reader;
    Json::Value root;

    std::ifstream input_file;
    input_file.open(input_path.c_str(), std::ios::binary);

    if (!reader.parse(input_file, root))
        throw std::runtime_error("Cannot parser json file.");

    if (!root["work_function"].isNull())
        _work_function = root["work_function"].asDouble();
    if (!root["fano"].isNull())
        _fano = root["fano"].asDouble();
    if (!root["transverse_diffusion"].isNull())
        _transverse_diffusion = root["transverse_diffusion"].asDouble();
    if (!root["longitudinal_diffusion"].isNull())
        _longitudinal_diffusion = root["longitudinal_diffusion"].asDouble();
    if (!root["drift_velocity"].isNull())
        _drift_velocity = root["drift_velocity"].asDouble();
    if (!root["z_plane"].isNull())
        _z_plane = root["z_plane"].asDouble();
    if (!root["trigger_threshold"].isNull())
        _trigger_threshold = root["trigger_threshold"].asDouble();
    if (!root["micromegas_size"].isNull())
        _micromegas_size = root["micromegas_size"].asDouble();
    if (root["micromegas_structure"].isNull())
        throw std::runtime_error("MicroMegas structure error.");
    else {
        int s = root["micromegas_structure"].size();
        for (int i = 0; i < s; i++) {
            double offset_y = (i + 0.5 - s/2.) * _micromegas_size;
            int k = root["micromegas_structure"][i].asInt();
            for (int j = 0; j < k; j++) {
                double offset_x = (j + 0.5 - k/2.) * _micromegas_size;
                _micromegas_info.push_back(make_tuple(offset_x, offset_y));
            }
        }
    }
}

int PrototypeConvert::get_electron_numbers(double energy)
{
    if (energy <= 0)
        return 0;
    int n = int(round(random.Gaus(energy / _work_function, TMath::Sqrt(_fano * energy / _work_function))));
    return n > 0 ? n : 0;
}

void PrototypeConvert::process(int i)
{
    chain->GetEntry(i);
    _gap_count = 0;
    _record_count = 0;
    _triggered = false;
    _trigger_energy = 0;
    _electron_info.clear();

    get_electron_info();
    auto result = get_trigger_info(_electron_info);
    if (get<0>(result) == -1) {
        _triggered = false;
        _trigger_energy = 0;
    } else {
        _triggered = true;
        _trigger_energy = (get<2>(result) - get<0>(result)) * _work_function;
    }
    _output_tree->Fill();
}

void PrototypeConvert::get_electron_info()
{
    if ((*xd).size() <= 0)
        return;
    for (auto i = 0; i < (*xd).size(); i++) {
        auto x = (*xd)[i];
        auto y = (*yd)[i];
        auto z = (*zd)[i];
        auto t = (*td)[i];
        auto e = (*energy)[i];
        double drift_z = _z_plane - z;
        if (drift_z < 0)
            continue;
        double drift_time = drift_z / _drift_velocity;//in us
        double r_diffusion_sigma = _transverse_diffusion * TMath::Sqrt(drift_z / 10) * 10;
        double z_diffusion_sigma = _longitudinal_diffusion * TMath::Sqrt(drift_z / 10) * 10;
        int electron_numbers = get_electron_numbers(e);
        for (auto j = 0; j < electron_numbers; j++) {
            double r_diffusion = random.Gaus(0, r_diffusion_sigma);
            double z_diffusion = random.Gaus(0, z_diffusion_sigma);
            double theta = random.Rndm() * 2 * TMath::Pi();
            double xx = x + r_diffusion * TMath::Sin(theta);
            double yy = y + r_diffusion * TMath::Cos(theta);
            double tt = t * 1e6 + drift_time + z_diffusion / _drift_velocity;
            if (out_of_border(xx, yy))
                _gap_count++;
            else {
                _record_count++;
                _electron_info.push_back(tt);
            }
        }
    }
}

bool PrototypeConvert::out_of_border(double x, double y)
{
    for (auto &u:_micromegas_info)
        if ((abs(x - get<0>(u)) < _micromegas_size/2.) && (abs(y - get<1>(u)) < _micromegas_size/2.))
            return false;
    return true;
}

tuple<int, int, int> PrototypeConvert::get_trigger_info(vector<double> &electron_info)
{
    sort(electron_info.begin(), electron_info.end(),
         [](double &l, double &r) { return l < r; });

    if (electron_info.size() <= 0)
        return make_tuple(-1, -1, -1);
    int start = 0;
    int electron_num_need = int(floor(_trigger_threshold / _work_function));
    bool find_it = false;
    while (start + electron_num_need < electron_info.size()) {
        if ((electron_info[start + electron_num_need] - electron_info[start]) <
            _trigger_bin_offset * _trigger_bin_length) {
            find_it = true;
            break;
        }
        start++;
    }
    if (find_it) {
        int trigger = start + electron_num_need;
        double end_time = electron_info[start + electron_num_need] +
                          (_trigger_bin_size - _trigger_bin_offset) * _trigger_bin_length;
        int end = int(find_if(electron_info.begin(), electron_info.end(),
                              [end_time](double l) { return l > end_time; }) - electron_info.begin());
        return make_tuple(start, trigger, end);
    }

    return make_tuple(-1, -1, -1);
}

void PrototypeConvert::final()
{
    process_all();
    _output_tree->Write();
    _output_file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-j", "--json", 1, false);
    parser.addArgument("-o", "--output", 1, true);
    parser.parse(argc, argv);
    PrototypeConvert *convert;
    if (parser.count("o"))
        convert = new PrototypeConvert(parser.retrieve<string>("input"), parser.retrieve<string>("json"),
                                       parser.retrieve<string>("output"));
    else
        convert = new PrototypeConvert(parser.retrieve<string>("input"), parser.retrieve<string>("json"),
                                       "output.root");
    convert->final();
}

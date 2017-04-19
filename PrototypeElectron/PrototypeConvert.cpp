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

class PrototypeConvert
        : public BaseConvert
{
public:
    PrototypeConvert(string &input_path, string &output_path, string &json_path);

    void load_parameter(string &input_path);

    void final();

    void process(int i);

private:
    int get_electron_numbers(double energy);

    void get_electron_info();

    TFile *_output_file;
    TTree *_output_tree;
    TRandom random;

    double work_function;
    double fano;
    double transverse_diffusion;
    double longitudinal_diffusion;
    double drift_velocity;
    double z_plane;
    double trigger_threshold;

    int _gap_count;
    int _record_count;
    bool _triggered;
    double _trigger_Energy;

    vector<double> _electron_info;
};

PrototypeConvert::PrototypeConvert(string &input_path, string &output_path, string &json_path)
        : BaseConvert(input_path, "mcTree"), random(), _electron_info()
{
    work_function = 21.9 / 1000;
    fano = 0.14;
    transverse_diffusion = 0.0101968; // cm1/2
    longitudinal_diffusion = 0.0139049; // cm1/2
    drift_velocity = 1.87139; // mm/us
    z_plane = 990; //mm
    trigger_threshold = 1200;//keV
    load_parameter(json_path);

    _output_file = new TFile(output_path.c_str(), "RECREATE");
    _output_tree = new TTree("DriftResult", "DriftResult");

    _output_tree->Branch("triggerEnergy", &_trigger_Energy, "triggerEnergy/D");
    _output_tree->Branch("triggered", &_triggered);
    _output_tree->Branch("gapCount", &_gap_count);
    _output_tree->Branch("recordCount", &_record_count);
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
        work_function = root["work_function"].asDouble();
    if (!root["fano"].isNull())
        fano = root["fano"].asDouble();
    if (!root["transverse_diffusion"].isNull())
        transverse_diffusion = root["transverse_diffusion"].asDouble();
    if (!root["longitudinal_diffusion"].isNull())
        longitudinal_diffusion = root["longitudinal_diffusion"].asDouble();
    if (!root["drift_velocity"].isNull())
        drift_velocity = root["drift_velocity"].asDouble();
    if (!root["z_plane"].isNull())
        z_plane = root["z_plane"].asDouble();
    if (!root["trigger_threshold"].isNull())
        trigger_threshold = root["trigger_threshold"].asDouble();
}

int PrototypeConvert::get_electron_numbers(double energy)
{
    if (energy <= 0)
        return 0;
    int n = int(round(random.Gaus(energy / work_function, TMath::Sqrt(fano * energy / work_function))));
    return n > 0 ? n : 0;
}

void PrototypeConvert::get_electron_info()
{
    _electron_info.clear();
    if ((*xd).size() <= 0)
        return;
    for (auto i = 0; i < (*xd).size(); i++) {
        auto x = (*xd)[i];
        auto y = (*yd)[i];
        auto z = (*zd)[i];
        auto t = (*td)[i];
        auto e = (*energy)[i];
        double drift_z = z_plane - z;
        if (drift_z < 0)
            continue;
        double drift_time = drift_z / drift_velocity;//in us
        double r_diffusion_sigma = transverse_diffusion * TMath::Sqrt(drift_z / 10) * 10;
        double z_diffusion_sigma = longitudinal_diffusion * TMath::Sqrt(drift_z / 10) * 10;
        int electron_numbers = get_electron_numbers(e);
        for (auto j = 0; j < electron_numbers; j++) {
            double r_diffusion = random.Gaus(0, r_diffusion_sigma);
            double z_diffusion = random.Gaus(0, z_diffusion_sigma);
            double theta = random.Rndm() * 2 * TMath::Pi();
            double xx = x + r_diffusion * TMath::Sin(theta);
            double yy = y + r_diffusion * TMath::Cos(theta);
            double tt = t * 1e6 + drift_time + z_diffusion / drift_velocity;
            int channel_id = plane.GetChannel(xx, yy);
            //int channel_id = 1;
            if (channel_id >= 0) {
                electron_info.push_back(make_tuple(channel_id, tt));
                record_count++;
            } else if (channel_id == -1)
                out_pixel_count++;
            else if (channel_id == -2)
                gap_count++;
        }
    }
}
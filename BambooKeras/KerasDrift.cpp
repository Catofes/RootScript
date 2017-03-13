//
// Created by herbertqiao on 2/28/17.
//

#include <iostream>
#include <TRandom.h>
#include <TMath.h>
#include <TFile.h>
#include <map>
#include "argparse.h"
#include "ReadoutStruct.hh"
#include "base_convert.h"
#include "ReadoutWave.hh"


using namespace std;

typedef vector<int> Wave;
typedef tuple<int, double> ElectronInfo; // channel_id,t

struct RawHits
{
    vector<double> xd, yd, zd, td, energy;
};

class ConvertEvent
{
public:
    ConvertEvent(const string &);

    ReadoutWave convert(const RawHits &input);

    int gap_count;
    int out_pixel_count;
    int record_count;
    bool cross_cathode;

private:
    ReadoutPlane plane;

    TRandom random;

    double work_function;
    double fano;
    double transverse_diffusion;
    double longitudinal_diffusion;
    double drift_velocity;
    double z_plane;
    double trigger_threshold;

    int get_electron_numbers(double energy);

    tuple<int, int, int> get_trigger_offset(vector<ElectronInfo> &electron_info);

};

ConvertEvent::ConvertEvent(const string &json_file)
        : plane(json_file), random(222)
{
    work_function = 21.9 / 1000;
    fano = 0.14;
    transverse_diffusion = 0.0101968; // cm1/2
    longitudinal_diffusion = 0.0139049; // cm1/2
    drift_velocity = 1.87139; // mm/us
    z_plane = 990; //mm
    trigger_threshold = 1200;//keV
}

ReadoutWave ConvertEvent::convert(const RawHits &input)
{
    gap_count = 0;
    out_pixel_count = 0;
    record_count = 0;
    cross_cathode = false;

    if (input.xd.size() <= 0)
        return ReadoutWave();
    bool in_top = input.zd[0] > 0;

    vector<ElectronInfo> electron_info;
    for (auto i = 0; i < input.xd.size(); i++) {
        auto x = input.xd[i];
        auto y = input.yd[i];
        auto z = input.zd[i];
        auto t = input.td[i];
        auto e = input.energy[i];
        if (z > 0 != in_top) {
            cross_cathode = true;
            return ReadoutWave();
        }
        double drift_z = z_plane - abs(z);
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
    sort(electron_info.begin(), electron_info.end(),
         [](ElectronInfo &l, ElectronInfo &r) { return get<1>(l) < get<1>(r); });

    auto trigger_info = get_trigger_offset(electron_info);
    if (get<0>(trigger_info) < 0) {
        return ReadoutWave();
    }
    ReadoutWave result;
    result.trigger_offset = get<1>(electron_info[get<1>(trigger_info)]);
    for (int i = get<0>(trigger_info); i < get<2>(trigger_info); i++) {
        int channel_id = get<0>(electron_info[i]);
        int t = int(floor((get<1>(electron_info[i]) - result.trigger_offset) / 0.2 + 256));
        if (result.detectors.find(channel_id) == result.detectors.end()) {
            Wave wave;
            wave.resize(512, 0);
            wave[t]++;
            result.detectors[channel_id] = wave;
        } else {
            result.detectors[channel_id][t]++;
        }
    }
    result.begin = get<0>(trigger_info);
    result.trigger = get<1>(trigger_info);
    result.end = get<2>(trigger_info);
    result.total_energy = (get<2>(trigger_info) - get<0>(trigger_info)) * work_function;
    return result;
}


tuple<int, int, int> ConvertEvent::get_trigger_offset(vector<ElectronInfo> &electron_info)
{
    if (electron_info.size() <= 0)
        return make_tuple(-1, -1, -1);
    int start = 0;
    int electron_num_need = int(floor(trigger_threshold / work_function));
    bool find_it = false;
    while (start + electron_num_need < electron_info.size()) {
        if ((get<1>(electron_info[start + electron_num_need]) - get<1>(electron_info[start])) < 256 * 0.2) {
            find_it = true;
            break;
        }
        start++;
    }
    if (find_it) {
        int trigger = start + electron_num_need;
        double end_time = get<1>(electron_info[start + electron_num_need]) + 256 * 0.2;
        int end = int(find_if(electron_info.begin(), electron_info.end(),
                              [end_time](ElectronInfo l) { return get<1>(l) > end_time; }) - electron_info.begin());
        return make_tuple(start, trigger, end);
    }

    return make_tuple(-1, -1, -1);
}

int ConvertEvent::get_electron_numbers(double energy)
{
    if (energy <= 0)
        return 0;
    int n = int(round(random.Gaus(energy / work_function, TMath::Sqrt(fano * energy / work_function))));
    return n > 0 ? n : 0;
}

class Raw2Electron
        : public BaseConvert
{
public:
    Raw2Electron(const string &input_path, const string &json_file, const string &output_path = "readout_wave.root");

    void final();

    void process(int i);

private:
    TFile *_output_file;
    TTree *_output_tree;
    ConvertEvent convert_event;
    ReadoutWave *_readout_wave;
    double _trigger_Energy;

    bool _triggered;
    bool _cross_cathode;
    int _gap_count;
    int _out_pixel_count;
    int _record_count;

};

Raw2Electron::Raw2Electron(const string &input_path, const string &json_file, const string &output_path)
        : BaseConvert(input_path), convert_event(json_file)
{
    _readout_wave = new ReadoutWave();
    _output_file = new TFile(output_path.c_str(), "RECREATE");
    _output_tree = new TTree("ReadoutWave", "ReadoutWave");
    _output_tree->Branch("runId", &runId, "runId/I");
    _output_tree->Branch("eventId", &eventId, "eventId/I");
    _output_tree->Branch("totalEnergy", &totalEnergy, "totalEnergy/D");
    _output_tree->Branch("uuid", &uuid);
    _output_tree->Branch("triggerEnergy", &_trigger_Energy, "triggerEnergy/D");
    _output_tree->Branch("readoutWave", &_readout_wave);
    _output_tree->Branch("triggered", &_triggered);
    _output_tree->Branch("gapCount", &_gap_count);
    _output_tree->Branch("outPixelCount", &_out_pixel_count);
    _output_tree->Branch("recordCount", &_record_count);
    _output_tree->Branch("crossCathode", &_cross_cathode);

}

void Raw2Electron::process(int i)
{
    if (i > 10000)
        return;
    if (i % 10 == 0)cout << i << endl;
    chain->GetEntry(i);
//    if (totalEnergy < 1200)
//        return;
    RawHits input;
    input.xd = *xd;
    input.yd = *zd;
    input.zd = *yd;
    input.td = *td;
    input.energy = *energy;
    auto event_result = convert_event.convert(input);
    //if (event_result.total_energy > 0) {
    _trigger_Energy = event_result.total_energy;
    _readout_wave = &event_result;
    _triggered = event_result.total_energy > 0;
    _gap_count = convert_event.gap_count;
    _out_pixel_count = convert_event.out_pixel_count;
    _record_count = convert_event.record_count;
    _cross_cathode = convert_event.cross_cathode;
    _output_tree->Fill();
    //}
}

void Raw2Electron::final()
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
    Raw2Electron *convert;
    if (parser.count("o"))
        convert = new Raw2Electron(parser.retrieve<string>("input"), parser.retrieve<string>("json"),
                                   parser.retrieve<string>("output"));
    else
        convert = new Raw2Electron(parser.retrieve<string>("input"), parser.retrieve<string>("json"));
    convert->final();
}
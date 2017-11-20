#include <argparse.h>
#include "ReadoutWave.hh"
#include "iostream"
#include "base_convert.h"
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <ReadoutStruct.hh>
#include <TMath.h>
#include <TLine.h>
#include <TText.h>
#include <TCanvas.h>
#include <sys/wait.h>

struct RawHits
{
    vector<double> xd, yd, zd, td, energy;
};

typedef vector<int> Wave;
typedef tuple<int, double> ElectronInfo; // channel_id,t

class T3DConvert
        : public BaseConvert
{
public:
    T3DConvert(const string &input_path, const string &json_file, const string &output_path);

    void process(int i);

    void final();

    void final(int offset, int limit, int pid);

    ReadoutWave raw2wave(const RawHits &input);

    void reconstruct(const ReadoutWave &input);

    tuple<double, double, double> get_center();

    void cluster();

    int get_electron_numbers(double energy);

    tuple<int, int, int> get_trigger_offset(vector<ElectronInfo> &electron_info);

    TFile *output_file;
    TTree *output_tree;
    double _trigger_Energy;
    bool _triggered;
    bool _cross_cathode;
    int _gap_count;
    int _out_pixel_count;
    int _record_count;
    vector<double> *_x = new vector<double>;
    vector<double> *_y = new vector<double>;
    vector<double> *_z = new vector<double>;
    vector<double> *_e = new vector<double>;

    vector<int> *_cluster_x = new vector<int>;
    vector<int> *_cluster_y = new vector<int>;
    vector<int> *_cluster_z = new vector<int>;
    vector<double> *_cluster_e = new vector<double>;
    string _type;
    double _cluster_max_energy;
    int _usage;


    int gap_count;
    int out_pixel_count;
    int record_count;
    bool cross_cathode;

    ReadoutPlane plane;

    TRandom random;

    double work_function;
    double fano;
    double transverse_diffusion;
    double longitudinal_diffusion;
    double drift_velocity;
    double z_plane;
    double trigger_threshold;

    const double cluster_x_size = 3;
    const double cluster_y_size = 3;
    const double cluster_z_size = 3;

};

T3DConvert::T3DConvert(const string &input_file_name, const string &json_file, const string &output_path)
        : BaseConvert(input_file_name, "PictureData"), plane(json_file), random(222)
{
    work_function = 21.9 / 1000;
    fano = 0.14;
    transverse_diffusion = 0.0101968; // cm1/2
    longitudinal_diffusion = 0.0139049; // cm1/2
    drift_velocity = 1.87139; // mm/us
    z_plane = 990; //mm
    trigger_threshold = 1200;//keV


    output_file = new TFile(output_path.c_str(), "RECREATE");
    output_tree = new TTree("3DPictureData", "Picture Data");

    output_tree->Branch("runId", &runId, "runId/I");
    output_tree->Branch("eventId", &eventId, "eventId/I");
    output_tree->Branch("nHits", &nHits, "nHits/I");
    output_tree->Branch("totalEnergy", &totalEnergy);
    output_tree->Branch("energy", &energy);
    output_tree->Branch("uuid", &uuid);

    output_tree->Branch("triggerEnergy", &_trigger_Energy);
    output_tree->Branch("triggered", &_triggered);
    output_tree->Branch("gapCount", &_gap_count);
    output_tree->Branch("outPixelCount", &_out_pixel_count);
    output_tree->Branch("recordCount", &_record_count);
    output_tree->Branch("crossCathode", &_cross_cathode);

    output_tree->Branch("x", &_x);
    output_tree->Branch("y", &_y);
    output_tree->Branch("z", &_z);
    output_tree->Branch("e", &_e);

    output_tree->Branch("cx", &_cluster_x);
    output_tree->Branch("cy", &_cluster_y);
    output_tree->Branch("cz", &_cluster_z);
    output_tree->Branch("ce", &_cluster_e);
    output_tree->Branch("cluster_max_energy", &_cluster_max_energy);
    output_tree->Branch("event_type", &_type);
    output_tree->Branch("usage", &_usage);

}

void T3DConvert::process(int i)
{
    chain->GetEntry(i);
    RawHits input;
    input.xd = *xd;
    input.yd = *zd;
    input.zd = *yd;
    input.td = *td;
    input.energy = *energy;
    auto result = raw2wave(input);
    reconstruct(result);
    cluster();
    output_tree->Fill();
}

void T3DConvert::final()
{
    process_all();
    output_tree->Write();
    output_file->Close();
}

void T3DConvert::final(int offset, int limit, int pid)
{
    cout << "Total " << limit << " Events." << endl;
    for (int i = offset; i < offset + limit; i++) {
        if (i % 1000 == 0)
            cout << "subprocess@" << pid << ":" << i << "/" << limit << endl;
        process(i);
    }
}

ReadoutWave T3DConvert::raw2wave(const RawHits &input)
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

void T3DConvert::reconstruct(const ReadoutWave &input)
{
    _x->clear();
    _y->clear();
    _z->clear();
    _e->clear();
    for (const auto &u:input.detectors) {
        int channel_id = u.first;
        double x = plane.GetX(channel_id);
        double y = plane.GetY(channel_id);
        for (int j = 0; j < u.second.size(); j++) {
            if (u.second[j] <= 0)
                continue;
            double e = u.second[j] * work_function;
            double z = j * 0.2 * drift_velocity;
            _x->push_back(x);
            _y->push_back(y);
            _z->push_back(z);
            _e->push_back(e);
        }
    }
}

tuple<double, double, double> T3DConvert::get_center()
{
    double total_energy = 0;
    double total_x = 0;
    double total_y = 0;
    double total_z = 0;
    for (int i = 0; i < _x->size(); i++) {
        total_energy += (*_e)[i];
        total_x += (*_x)[i] * (*_e)[i];
        total_y += (*_y)[i] * (*_e)[i];
        total_z += (*_z)[i] * (*_e)[i];
    }
    return make_tuple(total_x / total_energy, total_y / total_energy, total_z / total_energy);
};

void T3DConvert::cluster()
{
    _cluster_x->clear();
    _cluster_y->clear();
    _cluster_z->clear();
    _cluster_e->clear();
    map<tuple<int, int, int>, double> cluster_data;
    double center_x = 0, center_y, center_z = 0;
    tie(center_x, center_y, center_z) = get_center();
    for (int i = 0; i < (*_x).size(); i++) {
        auto cluster_x = (int) floor(((*_x)[i] - center_x) / cluster_x_size);
        auto cluster_y = (int) floor(((*_y)[i] - center_y) / cluster_y_size);
        auto cluster_z = (int) floor(((*_z)[i] - center_z) / cluster_z_size);
        auto cluster_id = make_tuple(cluster_x, cluster_y, cluster_z);
        if (cluster_data.find(cluster_id) != cluster_data.end())
            cluster_data[cluster_id] += (*_e)[i];
        else
            cluster_data.insert(make_pair(cluster_id, (*_e)[i]));
    }
    double max_energy = 0;
    for (auto &u:cluster_data) {
        if (u.second > max_energy)
            max_energy = u.second;
    }
    for (auto &u:cluster_data) {
        _cluster_x->push_back(get<0>(u.first));
        _cluster_y->push_back(get<1>(u.first));
        _cluster_z->push_back(get<2>(u.first));
        _cluster_e->push_back(u.second);
    }
    _cluster_max_energy = max_energy;
}

int T3DConvert::get_electron_numbers(double energy)
{
    if (energy <= 0)
        return 0;
    int n = int(round(random.Gaus(energy / work_function, TMath::Sqrt(fano * energy / work_function))));
    return n > 0 ? n : 0;
}

tuple<int, int, int> T3DConvert::get_trigger_offset(vector<ElectronInfo> &electron_info)
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

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-j", "--json", 1, false);
    parser.addArgument("-o", "--output", 1, true);
    parser.addArgument("-f", "--fork", 1, false);
    parser.parse(argc, argv);

    auto chain = new TChain("PictureData");
    chain->Add(parser.retrieve<string>("input").c_str());
    auto entries = chain->GetEntries();
    delete chain;

    int subprocess = stoi(parser.retrieve<string>("fork"));

    for (int i = 0; i < subprocess; i++) {
        auto child_pid = fork();

        if (child_pid == 0) {
            auto draw = new T3DConvert(parser.retrieve<string>("input"), parser.retrieve<string>("json"),
                                       parser.retrieve<string>("output") + to_string(i) + ".root");
            draw->final(i * (entries / subprocess), entries / subprocess, i);
        }
        if (child_pid < 0) {
            // Forking failed.
            perror("fork()");
            exit(EXIT_FAILURE);
        }
    }

    while (true) {
        int status;
        pid_t done = wait(&status);
        if (done == -1) {
            if (errno == ECHILD) break; // no more child processes
        } else {
            if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
                cerr << "pid " << done << " failed" << endl;
                exit(1);
            }
        }
    }

    //draw->final();
}
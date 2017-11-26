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

struct RawHits
{
    vector<double> xd, yd, zd, td, energy;
};

typedef vector<int> Wave;
typedef tuple<int, double> ElectronInfo; // channel_id,t

class DrawWave
        : public BaseConvert
{
public:
    DrawWave(const string &input_path, const string &json_file, const string &output_path);

    void process(int i);

    void final();

    pair<Wave, tuple<int, int, int>> convert(const RawHits &input);

    int get_electron_numbers(double energy);

    tuple<int, int, int> get_trigger_offset(vector<ElectronInfo> &electron_info);

    TFile *output_file;
    TH1F *output_th1f;
    TCanvas *canvas;
    ReadoutWave *_readout_wave = 0;
    double _trigger_Energy;
    bool _triggered;
    bool _cross_cathode;
    int _gap_count;
    int _out_pixel_count;
    int _record_count;
    vector<double> *_xzx = new vector<double>;
    vector<double> *_xzz = new vector<double>;
    vector<double> *_xze = new vector<double>;
    vector<double> *_yzy = new vector<double>;
    vector<double> *_yzz = new vector<double>;
    vector<double> *_yze = new vector<double>;

    vector<int> *_cluster_xzx = new vector<int>;
    vector<int> *_cluster_xzz = new vector<int>;
    vector<double> *_cluster_xze = new vector<double>;
    vector<int> *_cluster_yzy = new vector<int>;
    vector<int> *_cluster_yzz = new vector<int>;
    vector<double> *_cluster_yze = new vector<double>;

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

};

DrawWave::DrawWave(const string &input_file_name, const string &json_file, const string &output_path)
        : BaseConvert(input_file_name, "PictureData"), plane(json_file), random(222)
{
    work_function = 21.9 / 1000;
    fano = 0.14;
    transverse_diffusion = 0.0101968; // cm1/2
    longitudinal_diffusion = 0.0139049; // cm1/2
    drift_velocity = 1.87139; // mm/us
    z_plane = 990; //mm
    trigger_threshold = 1200;//keV

    chain->SetBranchAddress("triggerEnergy", &_trigger_Energy);
    chain->SetBranchAddress("readoutWave", &_readout_wave);
    chain->SetBranchAddress("triggered", &_triggered);
    chain->SetBranchAddress("gapCount", &_gap_count);
    chain->SetBranchAddress("outPixelCount", &_out_pixel_count);
    chain->SetBranchAddress("recordCount", &_record_count);
    chain->SetBranchAddress("crossCathode", &_cross_cathode);

    chain->SetBranchAddress("xzx", &_xzx);
    chain->SetBranchAddress("xzz", &_xzz);
    chain->SetBranchAddress("xze", &_xze);
    chain->SetBranchAddress("yzy", &_yzy);
    chain->SetBranchAddress("yzz", &_yzz);
    chain->SetBranchAddress("yze", &_yze);

    output_file = new TFile(output_path.c_str(), "RECREATE");
    canvas = new TCanvas("c1");
}

void DrawWave::process(int i)
{
    chain->GetEntry(i);
    RawHits input;
    input.xd = *xd;
    input.yd = *zd;
    input.zd = *yd;
    input.td = *td;
    input.energy = *energy;
    auto result = convert(input);
    auto wave = result.first;
    auto info = result.second;

    output_th1f = new TH1F("wave", "", wave.size(), 0, wave.size());
    output_th1f->GetXaxis()->SetTitle("Time bins");
    output_th1f->GetXaxis()->CenterTitle();
    output_th1f->GetYaxis()->SetTitle("Electron count");
    output_th1f->GetYaxis()->SetTitleOffset(1.2);
    output_th1f->GetXaxis()->CenterTitle();
    output_th1f->SetStats(false);

    for (auto k = 0; k < wave.size(); k++) {
        output_th1f->SetBinContent(k + 1, wave[k]);
    }

    output_th1f->Draw();

    auto y_max = output_th1f->GetMaximum();
    auto y_min = output_th1f->GetMinimum();

    TLine *line = new TLine(get<0>(info), y_max, get<0>(info), y_min);
    line->SetLineColor(kRed);
    line->Draw();

    TText *t = new TText(get<0>(info)+5, y_max, "Start Point");
    t->SetTextColor(kRed);
    t->SetTextSize(0.05);
    t->SetTextAngle(-90);
    t->Draw();

    line = new TLine(get<1>(info), y_max, get<1>(info), y_min);
    line->SetLineColor(kBlue);
    line->Draw();

    t = new TText(get<1>(info)+5, y_max, "Trigger Point");
    t->SetTextColor(kBlue);
    t->SetTextSize(0.05);
    t->SetTextAngle(-90);
    t->Draw();

    line = new TLine(get<2>(info), y_max, get<2>(info), y_min);
    line->SetLineColor(kMagenta);
    line->Draw();

    t = new TText(get<2>(info)-20, y_max, "End Point");
    t->SetTextColor(kMagenta);
    t->SetTextSize(0.05);
    t->SetTextAngle(-90);
    t->Draw();

    canvas->Update();
    canvas->Write();


}

void DrawWave::final()
{
    output_th1f->Write();
    output_file->Close();
}

pair<Wave, tuple<int, int, int>> DrawWave::convert(const RawHits &input)
{
    gap_count = 0;
    out_pixel_count = 0;
    record_count = 0;
    cross_cathode = false;

    if (input.xd.size() <= 0)
        return make_pair(Wave(), make_tuple(0, 0, 0));
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
            return make_pair(Wave(), make_tuple(0, 0, 0));
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

    Wave result;
    auto start_time = get<1>(electron_info[0]);
    auto trigger_point = int(floor((get<1>(electron_info[get<1>(trigger_info)]) - start_time) / 0.2));
    auto start_point = (trigger_point - 256) >= 0 ? trigger_point - 256 : 0;
    auto end_point = trigger_point + 256;
    auto max_t = int(floor((get<1>(electron_info[electron_info.size() - 1]) - start_time) / 0.2));
    if (max_t < end_point)
        max_t = end_point;
    result.resize(max_t + 1, 0);
    for (auto &hit:electron_info) {
        int bin = int(floor((get<1>(hit) - start_time) / 0.2));
        result[bin]++;
    }

    return make_pair(result, make_tuple(start_point, trigger_point, end_point));
}

int DrawWave::get_electron_numbers(double energy)
{
    if (energy <= 0)
        return 0;
    int n = int(round(random.Gaus(energy / work_function, TMath::Sqrt(fano * energy / work_function))));
    return n > 0 ? n : 0;
}

tuple<int, int, int> DrawWave::get_trigger_offset(vector<ElectronInfo> &electron_info)
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
    parser.addArgument("-e", "--entry", 1, false);
    parser.addArgument("-j", "--json", 1, false);
    parser.addArgument("-o", "--output", 1, true);
    parser.parse(argc, argv);
    auto draw = new DrawWave(parser.retrieve<string>("input"), parser.retrieve<string>("json"),
                             parser.retrieve<string>("output"));
    draw->process(atoi(parser.retrieve<string>("entry").c_str()));
    draw->final();
}
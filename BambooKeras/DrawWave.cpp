#include <argparse.h>
#include "ReadoutWave.hh"
#include "iostream"
#include "base_convert.h"
#include <TFile.h>
#include <TH1F.h>

class DrawWave
        : public BaseConvert
{
public:
    DrawWave(const string &input_path, const string &output_path);

    void process(int i);

    void final();

    TFile *output_file;
    TH1F *output_th1f;
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


};

DrawWave::DrawWave(const string &input_file_name, const string &output_path)
        : BaseConvert(input_file_name, "PictureData")
{
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
    output_th1f = new TH1F("wave", "Wave Form", 512, 0, 512);
    output_th1f->GetXaxis()->SetTitle("Time bins");
    output_th1f->GetYaxis()->SetTitle("Electron count");
}

void DrawWave::process(int i)
{
    chain->GetEntry(i);
    vector<int> results;
    for (auto i = 0; i < 512; i++)
        results.push_back(0);
    for (const auto &detector: _readout_wave->detectors) {
        auto wave = detector.second;
        for (auto i = 0; i < 512; i++) {
            results[i] += wave[i];
        }
    }
    for (auto i = 0; i < 512; i++) {
        output_th1f->SetBinContent(i + 1, results[i]);
    }
}

void DrawWave::final()
{
    output_th1f->Write();
    output_file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-e", "--entry", 1, false);
    parser.addArgument("-o", "--output", 1, true);
    parser.parse(argc, argv);
    auto draw = new DrawWave(parser.retrieve<string>("input"), parser.retrieve<string>("output"));
    draw->process(atoi(parser.retrieve<string>("entry").c_str()));
    draw->final();
}
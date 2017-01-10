#include "base_convert.h"
#include <TFile.h>

class AreaCutWithDrift
        : public BaseConvert
{
public:
    AreaCutWithDrift(const string &input_path, const double &rate, const string &output_path);

    void process(int i);

    void final();

private:
    TFile *out_file = 0;
    TH1F *h1 = 0;
    TH1F *h2 = 0;
    TH1F *h3 = 0;
    double _mean_free_path;
};

AreaCutWithDrift::AreaCutWithDrift(const string &input_path, const double &mean_free_path = 0,
                                   const string &output_path = "output.root")
        : BaseConvert(input_path)
{
    out_file = new TFile(output_path.c_str(), "RECREATE");
    h1 = new TH1F("energy_1", "Energy", 500, 0, 6000);
    h2 = new TH1F("energy_2", "Energy", 500, 0, 200);
    h3 = new TH1F("energy_3", "Energy", 500, 5000, 6000);
    _mean_free_path = mean_free_path;
}

void AreaCutWithDrift::process(int entry)
{
    chain->GetEntry(entry);
    if ((*primaryType)[0] != "Am241")
        return;
    double entry_energy = 0;
    for (int i = 0; i < xd->size(); i++) {
        double x = (*xd)[i];
        double y = (*yd)[i];
        double z = (*zd)[i];
        double e = (*energy)[i];
        string t = (*type)[i];
        if (t == "alpha" || t == "e-" || t == "gamma")
            if (-100 < x < 100 && -100 < y < 100) {
                if (_mean_free_path == 0)
                    entry_energy += e;
                else
                    entry_energy += e * exp(-(392.5 - z) / _mean_free_path);
            }
    }
    if (entry_energy == 0)
        return;
    h1->Fill(entry_energy);
    h2->Fill(entry_energy);
    h3->Fill(entry_energy);
}

void AreaCutWithDrift::final()
{
    process_all();
    out_file->cd();
    h1->Write();
    h2->Write();
    h3->Write();
    out_file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-r", "--rate", 1, false);
    parser.addArgument("-o", "--output", 1, false);
    parser.parse(argc, argv);
    TApplication *myapp = new TApplication("App", &argc, argv);
    AreaCutWithDrift convert(parser.retrieve<string>("input"), stod(parser.retrieve<string>("rate")),
                             parser.retrieve<string>("output"));
    cout << "Bind finished" << endl;
    convert.final();
    //myapp->Run();
}
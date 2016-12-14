#include <TFile.h>
#include <TChain.h>
#include <algorithm>
#include <iostream>
#include <TApplication.h>
#include "argparse.h"
#include <vector>
#include <tuple>
#include <TH1F.h>

using namespace std;

typedef tuple<double, double> hit;

bool compare(hit &left, hit &right)
{
    return get<1>(left) < get<1>(right);
}

bool draw(string input_path = "*.root", double time_cut = 1E-6, string input_tree_name = "mcTree")
{
    TChain input_chain(input_tree_name.c_str());
    input_chain.Add(input_path.c_str());

    vector<double> *e = 0, *td = 0;
    input_chain.SetBranchAddress("energy", &e);
    input_chain.SetBranchAddress("td", &td);

    vector<hit> hit_data;
    TH1F *h = new TH1F("energy_spectrum", "Energy Spectrum", 500, 0, 700);

    int n_entries = input_chain.GetEntries();
    cout << "Load " << n_entries << " entries." << endl;

    for (int i = 0; i < n_entries; i++) {
        input_chain.GetEntry(i);
        if (i % 1000 == 0)
            cout << i << "/" << n_entries << endl;
        hit_data.clear();
        hit_data.reserve(e->size());
        for (int j = 0; j < e->size(); j++) {
            hit_data.push_back(make_tuple((*e)[j], (*td)[j]));
        }
        sort(hit_data.begin(), hit_data.end(), compare);

        if (hit_data.size() == 0)
            continue;
        double total_energy = 0;
        double temp_time = get<1>(hit_data[0]);
        for (int j = 0; j < hit_data.size(); j++) {
            if ((get<1>(hit_data[j]) - temp_time) < time_cut)
                total_energy += get<0>(hit_data[j]);
            else {
                h->Fill(total_energy);
                temp_time = get<1>(hit_data[j]);
                total_energy = get<0>(hit_data[j]);
            }
        }
        h->Fill(total_energy);
    }
    h->Draw();
    return true;
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-t", "--time", 1);
    parser.parse(argc, argv);
    TApplication *myapp = new TApplication("App", &argc, argv);
    if (parser.count("time"))
        draw(parser.retrieve<string>("input"), parser.retrieve<float>("t"));
    else
        draw(parser.retrieve<string>("input"));
    myapp->Run();
}
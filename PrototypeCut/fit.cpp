//
// Created by herbertqiao on 3/2/17.
//

#include <iostream>
#include <TChain.h>
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include "argparse.h"

using namespace std;

class Fit
{
public:
    Fit(const string &input_file, const string &tree_name = "mcTree");

    void create_hist(const string &name = "totalEnergy");

    void load(const string &input_branch, const string &filter = "");

    void process(const string &output_file);

    int hist_bins = 500;
    double hist_start = 0;
    double hist_end = 70;
    double fit_start = 0;
    double fit_end = 70;
private:
    TChain *chain;
    TH1F *hist;
    TFile *file;
};

Fit::Fit(const string &input_file, const string &tree_name)
{
    chain = new TChain(tree_name.c_str());
    chain->Add(input_file.c_str());
}

void Fit::create_hist(const string &name)
{
    hist = new TH1F("hist", name.c_str(), hist_bins, hist_start, hist_end);
}

void Fit::load(const string &input_branch, const string &filter)
{
    string exp(input_branch);
    exp += ">>hist";
    chain->Draw(exp.c_str(), filter.c_str());
}

void Fit::process(const string &output_file)
{
    file = new TFile(output_file.c_str(), "RECREATE");
    TF1 *f1 = new TF1("f1", "gaus", fit_start, fit_end);
    hist->Fit("f1", "R");
    file->cd();
    hist->Write();
    file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.addArgument("-t", "--tree", 1, true);
    parser.addArgument("-b", "--branch", 1, false);
    parser.addArgument("-f", "--filter", 1, true);
    parser.addArgument("--hits_name", 1, true);
    parser.addArgument("--hits_bins", 1, true);
    parser.addArgument("--hist-start", 1, true);
    parser.addArgument("--hist-end", 1, true);
    parser.addArgument("--fit-start", 1, true);
    parser.addArgument("--fit-end", 1, true);
    parser.addArgument("-o", "--output", 1, false);

    parser.parse(argc, argv);
    Fit *fit;
    if (parser.count("tree"))
        fit = new Fit(parser.retrieve<string>("input"), parser.retrieve<string>("tree"));
    else
        fit = new Fit(parser.retrieve<string>("input"));
    if (parser.count("hist-start"))
        fit->hist_start = stod(parser.retrieve<string>("hist-start"));
    if (parser.count("hist-end"))
        fit->hist_start = stod(parser.retrieve<string>("hist-end"));
    if (parser.count("fit-start"))
        fit->hist_start = stod(parser.retrieve<string>("fit-start"));
    if (parser.count("fit-end"))
        fit->hist_start = stod(parser.retrieve<string>("fit-end"));
    if (parser.count("hits_bins"))
        fit->hist_bins = stoi(parser.retrieve<string>("hits_bins"));
    if (parser.count("hits_name"))
        fit->create_hist(parser.retrieve<string>("hist_name"));
    else
        fit->create_hist();
    if (parser.count("filter"))
        fit->load(parser.retrieve<string>("branch"), parser.retrieve<string>("filter"));
    else
        fit->load(parser.retrieve<string>("branch"));
    fit->process(parser.retrieve<string>("output"));
}
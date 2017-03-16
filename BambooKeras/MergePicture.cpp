//
// Created by herbertqiao on 3/16/17.
//

#include <iostream>
#include <TChain.h>
#include <TFile.h>
#include <argparse.h>

using namespace std;

void convert(const vector<string> &input_file_name, const string &output_file_name)
{
    TChain *chain = new TChain("MLPicture");
    for (auto &u:input_file_name)
        chain->Add(u.c_str());
    TFile *file = new TFile(output_file_name.c_str(), "RECREATE");
    chain->Write();
    file->Close();
}

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-o", "--output", 1, false);
    parser.addArgument("--input", "+", false);
    parser.parse(argc, argv);
    convert(parser.retrieve<vector<string>>("input"), parser.retrieve<string>("output"));
}
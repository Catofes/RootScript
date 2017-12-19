//
// Created by herbertqiao on 12/5/17.
//

#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TFitResult.h>
#include <argparse.h>
#include <string>

using namespace std;

int fit(string input_path)
{
    TFile *f = new TFile(input_path.c_str(), "READ");
    TH1F *h = (TH1F *) f->Get("time1");
    auto result = h->Fit("gaus", "RS", "+", 275, 325).Get()->Parameters();

    cout << "Result: " << input_path << "\t" << to_string(result[2]) << endl;
}

int main(int argc, char *argv[])
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1, false);
    parser.parse(argc, argv);
    fit(parser.retrieve<string>("input"));
}
//
// Created by herbertqiao on 11/13/16.
//

#ifndef ROOTSCRIPT_HIT2CHAIN_H
#define ROOTSCRIPT_HIT2CHAIN_H

#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <string.h>
#include <TFile.h>
#include <vector>
#include <stdlib.h>

using namespace std;

class Hit2Data
{
public:
    Hit2Data(const string &input_file, const string &output_file, int type);

    void Convert();

private:
    void Bind();

    TChain _input_chain;
    TFile _output_file;
    TTree _output_tree;

    vector<double> *_x = 0;
    vector<double> *_y = 0;
    vector<double> *_z = 0;
    vector<double> *_e = 0;

    double _total_energy;
    int _type;
};

#endif //ROOTSCRIPT_HIT2CHAIN_H

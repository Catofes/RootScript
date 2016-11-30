//
// Created by herbertqiao on 11/13/16.
//

#include <stdlib.h>
#include "hit2data.h"


Hit2Data::Hit2Data(const string &input_file, const string &output_file, int type = 0)
        : _input_chain("mcTree"),
          _output_file(output_file.c_str(), "RECREATE"),
          _output_tree("MLData", "MachineLearning Data"),
          _type(type)
{
    _input_chain.Add(input_file.c_str());
    Bind();
}

void Hit2Data::Bind()
{
    _input_chain.SetBranchAddress("xd", &_x);
    _input_chain.SetBranchAddress("yd", &_y);
    _input_chain.SetBranchAddress("zd", &_z);
    _input_chain.SetBranchAddress("energy", &_e);
    _input_chain.SetBranchAddress("totalEnergy", &_total_energy);

    _input_chain.GetEntry(0);

    _output_tree.Branch("x", &_x);
    _output_tree.Branch("y", &_y);
    _output_tree.Branch("z", &_z);
    _output_tree.Branch("e", &_e);
    _output_tree.Branch("type", &_type);
}

void Hit2Data::Convert()
{
    long n_entries = _input_chain.GetEntries();
    for (int i = 0; i < n_entries; i++) {
        if (i % 10000 == 0)
            cout << "Process" << i << "/" << n_entries << endl;
        _input_chain.GetEntry(i);
        if (_total_energy > 2457.83 + 125.249 / 2. || _total_energy < 2457.83 - 125.249 / 2.)
            continue;
        _output_tree.Fill();
    }
    _output_tree.Write();
    _output_file.Close();
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        throw runtime_error("Use Hit2Data input_file, output_file, type.");
    }
    Hit2Data hit2data(argv[1], argv[2], atoi(argv[3]));
    hit2data.Convert();
}
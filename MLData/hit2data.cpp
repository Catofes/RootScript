//
// Created by herbertqiao on 11/13/16.
//


#include "hit2data.h"


Hit2Data::Hit2Data(const string &input_file, const string &output_file, int type = 0)
        : _input_chain("mcTree"),
          _output_file(output_file.c_str(), "w"),
          _output_tree("MLData", "MachineLearning Data"),
          _type(_type)
{
    _input_chain.Add(input_file.c_str());
    Bind();
}

void Hit2Data::Bind()
{
    _input_chain.SetBranchAddress("xd", &_x);
    _input_chain.SetBranchAddress("yd", &_y);
    _input_chain.SetBranchAddress("zd", &_z);
    _input_chain.SetBranchAddress("ed", &_e);

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
        _input_chain.GetEntry(i);
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
    Hit2Data hit2data(argv[1], argv[2], int(argv[3]));
    hit2data.Convert();
}
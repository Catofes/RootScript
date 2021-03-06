//
// Created by herbertqiao on 7/6/16.
//

#ifndef ROOTSCRIPT_SMEARMC_HH
#define ROOTSCRIPT_SMEARMC_HH

#include <iostream>
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>

class SmearMC
{
public:
    SmearMC(std::string input_path, std::string output_path);

    void smear();

private:
    double _total_energy = 0;
    double _cha_mca_energy = 0;
    double _cha_energy = 0;
    double _scale = 0.0667788;

    TChain *_input_chain;

    TFile *_output_file;

    TTree *_output_tree;

    int _simulate_num = 1000000;

};

#endif //ROOTSCRIPT_SMEARMC_HH

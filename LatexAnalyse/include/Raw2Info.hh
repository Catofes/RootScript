//
// Created by herbertqiao on 6/28/16.
//

#ifndef ROOTSCRIPT_RAW2INFO_HH
#define ROOTSCRIPT_RAW2INFO_HH

#include <iostream>
#include <TFile.h>
#include <TTree.h>

using namespace std;

class Hit
{
public:
    Hit()
    { }

public:
    int _trackId;
    std::vector<int> _parentId;
    std::vector<std::string > _type;
    std::vector<std::string > _parentType;
    std::vector<std::string > _creatorProcess;
    std::vector<std::string > _depositionProcess;
    std::vector<std::string > _volume;
    std::vector<double> _xd;
    std::vector<double> _yd;
    std::vector<double> _zd;
    std::vector<double> _td;
    std::vector<double> _energy;
};


class Raw2Info
{
public:
    Raw2Info(string input);

    void Process(int entries);

private:
    TFile *file;
    TTree *tree;
    int entries;
    vector<double> energy;
    vector<double> xd;
    vector<double> yd;
    vector<double> zd;
    vector<double> td;
    vector<string> volume;
    vector<string> primaryType;
    vector<double> *_energy = &energy;
    vector<double> *_xd = &xd;
    vector<double> *_yd = &yd;
    vector<double> *_zd = &zd;
    vector<double> *_td = &td;
    vector<string> *_volume = &volume;
    vector<string> *_primaryType = &primaryType;
};


#endif //ROOTSCRIPT_RAW2INFO_HH

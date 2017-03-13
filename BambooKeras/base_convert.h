//
// Created by herbertqiao on 1/13/17.
//

#ifndef ROOTSCRIPT_LOAD_FILE_H
#define ROOTSCRIPT_LOAD_FILE_H

#include <TChain.h>
#include <iostream>
#include "argparse.h"
#include <vector>
#include <TH1F.h>
#include <TApplication.h>

using namespace std;

class BaseConvert
{
public:

    BaseConvert(const string &input_path);

    virtual void process(int i)
    {};

    virtual void final()
    {};

    virtual void process_all();

protected:
    TChain *chain;
    Int_t runId;
    Int_t eventId;
    Int_t nHits;
    vector<int> *trackId;
    vector<int> *parentId;
    vector<string> *type;
    vector<string> *parentType;
    vector<string> *creatorProcess;
    vector<string> *depositionProcess;
    vector<string> *volume;
    Double_t totalEnergy;
    vector<double> *xd;
    vector<double> *yd;
    vector<double> *zd;
    vector<double> *td;
    vector<double> *energy;
    Int_t nPrimaries;
    vector<string> *primaryType;
    vector<int> *primaryId;
    vector<double> *primaryEnergy;
    vector<double> *primaryPx;
    vector<double> *primaryPy;
    vector<double> *primaryPz;
    vector<double> *primaryX;
    vector<double> *primaryY;
    vector<double> *primaryZ;
    string uuid;
};

BaseConvert::BaseConvert(const string &input_path = "*.root")
{
    chain = new TChain("mcTree");
    chain->Add(input_path.c_str());

    trackId = 0;
    parentId = 0;
    type = 0;
    parentType = 0;
    creatorProcess = 0;
    depositionProcess = 0;
    volume = 0;
    xd = 0;
    yd = 0;
    zd = 0;
    td = 0;
    energy = 0;
    primaryType = 0;
    primaryId = 0;
    primaryEnergy = 0;
    primaryPx = 0;
    primaryPy = 0;
    primaryPz = 0;
    primaryX = 0;
    primaryY = 0;
    primaryZ = 0;

    chain->SetBranchAddress("runId", &runId);
    chain->SetBranchAddress("eventId", &eventId);
    chain->SetBranchAddress("nHits", &nHits);
    chain->SetBranchAddress("trackId", &trackId);
    chain->SetBranchAddress("parentId", &parentId);
    chain->SetBranchAddress("type", &type);
    chain->SetBranchAddress("parentType", &parentType);
    chain->SetBranchAddress("creatorProcess", &creatorProcess);
    chain->SetBranchAddress("depositionProcess", &depositionProcess);
    chain->SetBranchAddress("volume", &volume);
    chain->SetBranchAddress("totalEnergy", &totalEnergy);
    chain->SetBranchAddress("xd", &xd);
    chain->SetBranchAddress("yd", &yd);
    chain->SetBranchAddress("zd", &zd);
    chain->SetBranchAddress("td", &td);
    chain->SetBranchAddress("energy", &energy);
    chain->SetBranchAddress("nPrimaries", &nPrimaries);
    chain->SetBranchAddress("primaryType", &primaryType);
    chain->SetBranchAddress("primaryId", &primaryId);
    chain->SetBranchAddress("primaryEnergy", &primaryEnergy);
    chain->SetBranchAddress("primaryPx", &primaryPx);
    chain->SetBranchAddress("primaryPy", &primaryPy);
    chain->SetBranchAddress("primaryPz", &primaryPz);
    chain->SetBranchAddress("primaryX", &primaryX);
    chain->SetBranchAddress("primaryY", &primaryY);
    chain->SetBranchAddress("primaryZ", &primaryZ);
    chain->SetBranchAddress("uuid", &uuid);
}

void BaseConvert::process_all()
{
    int entries = chain->GetEntries();
    cout << "Total " << entries << " Events." << endl;
    for (int i = 0; i < entries; i++) {
        if (i % 1000 == 0)
            cout << i << "/" << entries << endl;
        process(i);
    }
}

#endif //ROOTSCRIPT_LOAD_FILE_H

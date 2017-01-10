//
// Created by herbertqiao on 12/21/16.
//

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <TChain.h>
#include <TMath.h>

using namespace std;

class ReadData
{
public:
    ReadData(string path);

    void print_distance();

    void save_distance();

private:
    double get_a_distance(int entry);

    TChain *chain;
    vector<double> *x;
    vector<double> *y;
    vector<double> *z;
    vector<double> *weight;
};

//This is the construct function. It read root file from path and bind the branch address.
ReadData::ReadData(string path)
{
    //Load the tree.
    chain = new TChain("Data");
    chain->Add(path.c_str());
    //Prepare the vector.
    x = 0, y = 0, z = 0, weight = 0;
    //Bind the branch.
    chain->SetBranchAddress("x", &x);
    chain->SetBranchAddress("y", &y);
    chain->SetBranchAddress("z", &z);
    chain->SetBranchAddress("weight", &weight);
}

//This function return the mass center's distance for a entry.
double ReadData::get_a_distance(int entry)
{
    int particle_number = x->size();
    if (particle_number <= 0)
        return 0;
    double total_weight = 0;
    double total_weight_multiply_distance = 0;
    for (int i = 0; i < particle_number; i++) {
        double distance = TMath::Sqrt((*x)[i] * (*x)[i] + (*y)[i] * (*y)[i] + (*z)[i] * (*z)[i]);
        total_weight += (*weight)[i];
        total_weight_multiply_distance += (*weight)[i] * distance;
    }
    return total_weight_multiply_distance / total_weight;
}

//This function print all entry's distance
void ReadData::print_distance()
{
    for (int i = 0; i < 100; i++)
        cout << "Center of entry " << i << " :" << get_a_distance(i) << endl;
}

//You need file 
void ReadData::save_distance()
{

}

int main()
{
    ReadData *read_data = new ReadData("input_data.root");
    read_data->print_distance();
}
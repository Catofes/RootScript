//
// Created by herbertqiao on 12/13/16.
//

#ifndef ROOTSCRIPT_G4DATACONVERT_H
#define ROOTSCRIPT_G4DATACONVERT_H

#include <iostream>
#include <string>
#include <TChain.h>
#include <vector>
#include <tuple>

using namespace std;

//Hit info with x, y, z, e at(mm, mm, mm, eV)
typedef tuple<double, double, double, double> hit;


class G4DataConvert
{
public:
    G4DataConvert(const string &path);

    ~G4DataConvert();

    int get_total_entries();

    vector<hit> get_hits(int entry);

    vector<double> *get_hit_x(int entry);

    vector<double> *get_hit_y(int entry);

    vector<double> *get_hit_z(int entry);

    vector<double> *get_hit_e(int entry);

private:
    TChain *chain = 0;
    vector<double> *x;
    vector<double> *y;
    vector<double> *z;
    vector<double> *e;

};


#endif //ROOTSCRIPT_G4DATACONVERT_H

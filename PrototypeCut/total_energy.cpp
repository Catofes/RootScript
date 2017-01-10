//
// Created by herbertqiao on 1/10/17.
//

#include <TChain.h>
#include <iostream>
#include "argparse.h"

using namespace std;

void draw(string input_path = "*.root")
{
    TChain *chain = new TChain("mcTree");
    chain->Add(input_path.c_str());

}
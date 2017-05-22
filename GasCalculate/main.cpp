#include <iostream>
#include "argparse.h"
#include <MediumMagboltz.hh>
#include <FundamentalConstants.hh>
#include <string.h>

using namespace std;
using namespace Garfield;

int main(int argc, char * argv[]) {
	ArgumentParser parser;
	parser.addArgument("-p", "--pressure", 1, false);
	parser.parse(argc, argv);
    const double pressure = stoi(parser.retrieve<string>("pressure")) / 1.01325 * AtmosphericPressure;
	const double temperature = 293.15;
    cout<<"Pressure: "<<pressure<<endl;
    // Setup the gas.
    MediumMagboltz* gas = new MediumMagboltz();
    gas->SetTemperature(temperature);
    gas->SetPressure(pressure);
    gas->SetComposition("Ar", 95., "ISO", 5);

    // Set the field range to be covered by the gas table.
    const int nFields = 1;
    const double emin = 200.;
    const double emax = 200.;
    // Flag to request logarithmic spacing.
    const bool useLog = false;
    gas->SetFieldGrid(emin, emax, nFields, useLog);

    const int ncoll = 10;
    // Switch on debugging to print the Magboltz output.
    gas->EnableDebugging();
    // Run Magboltz to generate the gas table.
    gas->GenerateGasTable(ncoll);
    gas->DisableDebugging();
    // Save the table.
    gas->WriteGasFile("xe_tma.gas");

}

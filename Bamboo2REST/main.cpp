#include "G4DataConvert.h"
#include "argparse.h"

int main(int argc, char **argv)
{
    ArgumentParser parser;
    parser.addArgument("-i", "--input", 1);
    parser.parse(argc, argv);
    //Load the data.
    //Parameter is roots' file path. Support "/foo/bar/*.root"
    G4DataConvert convert(parser.retrieve<string>("input"));

    //Get totally entries
    cout << "Total Entries: " << convert.get_total_entries() << endl;

    //Get the first entry
    vector<hit> fist_entry = convert.get_hits(0);
    cout << "First entry have " << fist_entry.size() << " hits." << endl;

    //Get the second entry's energy info. May faster than get_hits();
    vector<double> *second_entry_energy = convert.get_hit_e(1);
    double total_energy = 0;
    for (auto energy:*second_entry_energy)
        total_energy += energy;
    cout << "Second entry deposited " << total_energy << " keV." << endl;
}
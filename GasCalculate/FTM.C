#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TFile.h>
#include <TRandom3.h>

#include "ComponentAnsys123.hh"
#include "ViewField.hh"
#include "TrackHeed.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewDrift.hh"
#include "DriftLineRKF.hh"
#include "ViewCell.hh"
#include "ViewGeometry.hh"
#include "ViewSignal.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentAnalyticField.hh"
#include "ViewFEMesh.hh"

#include "argparse.h"

using namespace Garfield;
using namespace std;
namespace var
{
    double tau;
    string input_file_path;
    string field_path;
    string field_num;
    string gas_folder;
    string gas_name;
    string output_path;
    string layer_num;
}


double transfer(double t)
{
    const double tau = var::tau;
    return (t / tau) * (t / tau) * exp(1 - 2 * t / tau);
}

int shiyan()
{
    int nevents;
    int gapofele;
    double momentum;
    std::string particlename;
    double track_x;
    double dx;
    double dy;
    double dz;

    string _a, _b, _c, _d, _e;

    ifstream inp(var::input_file_path.c_str());
    inp >> nevents >> gapofele >> momentum >> particlename >> track_x >> dx >> dy >> dz;
    std::cout << "number of events:  " << nevents << "\ngap of events(ns):  " << gapofele
              << "\nenergy of particle(ev):  " << momentum << "\nparticle name:  " << particlename
              << "\nposition in X-axis(cm): " << track_x << "\ndx:  " << dx << "\ndy:  " << dy << "\ndz:  " << dz
              << std::endl;
    inp.close();

    int drawbin = nevents * gapofele;
    TH1D *h1 = new TH1D("h1", "", drawbin, 0, drawbin);
    TH1D *h2 = new TH1D("h2", "", drawbin, 0, drawbin);
    TH1D *h3 = new TH1D("h3", "", drawbin, 0, drawbin);
    TH1D *h4 = new TH1D("h4", "", drawbin, 0, drawbin);
    TH1D *h_iniele = new TH1D("h_iniele", "", 3000, 0, 3000);
    h_iniele->GetXaxis()->SetTitle("initial_electron");

    for (int pnum = 0; pnum < nevents; pnum++) {

        //注意修改电场文件路径
        ComponentAnsys123 *fm = new ComponentAnsys123();
        _a = var::field_path + "/" + var::field_num + "/ELIST.lis";
        _b = var::field_path + "/" + var::field_num + "/NLIST.lis";
        _c = var::field_path + "/" + var::field_num + "/MPLIST.lis";
        _d = var::field_path + "/" + var::field_num + "/PRNSOL.lis";
        _e = var::field_path + "/" + var::field_num + "/weight1.lis";
//        fm->Initialise("./field_1100/ELIST.lis", "./field_1100/NLIST.lis", "./field_1100/MPLIST.lis",
//                       "./field_1100/PRNSOL.lis", "mm");
//        fm->SetWeightingField("./field_1100/weight1.lis", "readout1");
        fm->Initialise(_a, _b, _c, _d, "mm");
        fm->SetWeightingField(_e, "readout1");

        fm->EnableMirrorPeriodicityX();
        fm->EnableMirrorPeriodicityY();

        //注意修改气体组分
        MediumMagboltz *gas = new MediumMagboltz();
        _a = var::gas_folder + "/" + var::gas_name + ".gas";
        gas->LoadGasFile(_a);


        gas->SetMaxElectronEnergy(200.);
        gas->Initialise(true);

        const double rPenning = 0.57;
        const double lambdaPenning = 0.;
        gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
        _a = var::gas_folder + "/IonMobility_Ar+_Ar.txt";
        gas->LoadIonMobility(_a);


        const int nMaterials = fm->GetNumberOfMaterials();
        for (int i = 0; i < nMaterials; ++i) {
            const double eps = fm->GetPermittivity(i);
            if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);
        }
        fm->PrintMaterials();


        const double pitch = 0.014;
        const double kapton = 50.e-4;
        const double metal = 5.e-4;
        const double outdia = 70.e-4;
        const double middia = 50.e-4;


        Sensor *sensor = new Sensor();
        sensor->SetTransferFunction(transfer);
        sensor->AddComponent(fm);
        sensor->AddElectrode(fm, "readout1");
        sensor->SetArea(0 * pitch, -pitch, 0.0,
                        3 * pitch, pitch, 0.09);


        AvalancheMicroscopic *aval = new AvalancheMicroscopic();
        aval->SetSensor(sensor);
        aval->EnableSignalCalculation();

        AvalancheMC *drift = new AvalancheMC();
        drift->SetSensor(sensor);
        drift->SetDistanceSteps(2.e-4);


        const double tmin = 0.;
        const double tmax = gapofele;
        const double tstep = 1;
        const int nTimeBins = int((tmax - tmin) / tstep);
        sensor->SetTimeWindow(0, tstep, nTimeBins);


        TrackHeed *track = new TrackHeed();
        track->SetParticle(particlename);
        track->SetMomentum(momentum);
        track->SetSensor(sensor);
        track->EnableElectricField();

        double xcls, ycls, zcls, tcls, e, extra;
        xcls = ycls = zcls = tcls = e = extra = -999.;

        int n = 0;

        //    double track_x =0.02;
        double track_y = 0.0;
        double track_z = 0.0788;

        double track_dx = dx;
        double track_dy = dy;
        double track_dz = dz;


        int numberofelectron = 0;
        double xele, yele, zele, tele, eele, dxele, dyele, dzele;

        int ionelectron = 0;
        track->NewTrack(track_x, track_y, track_z, tmin, track_dx, track_dy, track_dz);

        while (track->GetCluster(xcls, ycls, zcls, tcls, n, e, extra)) {
            std::cout << "number of event: " << pnum + 1 << "      clustersize=" << n << std::endl;
            for (int j = 1; j <= n; j++) {
                track->GetElectron(j - 1, xele, yele, zele, tele, eele, dxele, dyele, dzele);
                aval->AvalancheElectron(xele, yele, zele, tele, eele, dxele, dyele, dzele);


                unsigned int np = aval->GetNumberOfElectronEndpoints();
                double xe1, ye1, ze1, te1, e1;
                double xe2, ye2, ze2, te2, e2;
                double xi1, yi1, zi1, ti1;
                double xi2, yi2, zi2, ti2;

                int status1ele;
                for (int i = np; i--;) {
                    aval->GetElectronEndpoint(i, xe1, ye1, ze1, te1, e1, xe2, ye2, ze2, te2, e2, status1ele);
                }
                numberofelectron++;
                ionelectron++;
            }
        }
        std::cout << "number of initial electrons in this event: " << ionelectron << std::endl;
        h_iniele->Fill(ionelectron);
        sensor->ConvoluteSignal();
        for (int i = 0; i < gapofele; i++) {
            double a = sensor->GetSignal("readout1", i + 1);
            h1->SetBinContent(pnum * gapofele + i + 1, a);
        }

        delete fm;
        delete gas;
        delete sensor;
        delete aval;
        delete drift;
        delete track;

    }

    _a = var::output_path + "/FTM_L" + var::layer_num + "_" + var::gas_name + "_" + particlename
         + "_F" + var::field_num + "_T" + to_string(var::tau) + ".root";
    TFile *outfile1 = new TFile(_a.c_str(), "RECREATE");
    h1->Write();
    h_iniele->Write();
    outfile1->Close();

    std::cout << "Congratulations! ^-^ The data has been all saved! Please exit the program!" << std::endl;
}


int main(int argc, char *argv[])
{
//	TApplication theApp("App", &argc, argv);
    ArgumentParser parser;
    parser.addArgument("-t", "--tau", 1, false);
    parser.addArgument("-", "--gas_name", 1, false);
    parser.addArgument("-", "--gas_path", 1, false);
    parser.addArgument("-", "--input", 1, false);
    parser.addArgument("-", "--field_path", 1, false);
    parser.addArgument("-", "--field_num", 1, false);
    parser.addArgument("-", "--layer_num", 1, false);
    parser.addArgument("-", "--output", 1, false);

    parser.parse(argc, argv);
    var::tau = std::stod(parser.retrieve<string>("tau"));
    var::gas_name = parser.retrieve<string>("gas_name");
    var::gas_folder = parser.retrieve<string>("gas_folder");
    var::input_file_path = parser.retrieve<string>("input");
    var::field_path = parser.retrieve<string>("field_path");
    var::field_num = parser.retrieve<string>("field_num");
    var::layer_num = parser.retrieve<string>("layer_num");
    var::output_path = parser.retrieve<string>("output");


    shiyan();
    //	theApp.Run();
    return 0;
}

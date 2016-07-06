/*
  this kumac smears the name-switched MC ntuple
  - gesim Geant4 simulation
  - name-switching
  - smearing MC energy 
     using formula sigma = 0.512+0.0002632*E [keV]
  - compare with data

  this kumac will smear the Cha_Energy and changes the histograms
  MC always assumes 1Bq
*/

void smear_mc_on_bluster()
{

    const Int_t nisotopes = 7;
    Char_t *isotopename[nisotopes];
/*    
    isotopename[0]="rn222";
    isotopename[1]="rn220";
*/
    isotopename[0] = "u238";
    isotopename[1] = "th232";
    isotopename[2] = "k40";
    isotopename[3] = "co60";
    isotopename[4] = "u235";
    isotopename[5] = "cs137";
    isotopename[6] = "pb210";
/*
    isotopename[7]="ag110";
    isotopename[8]="ag110m";

    isotopename[0]="co56";
    isotopename[1]="co58";
    isotopename[2]="mn54";
    isotopename[3]="be7";
    isotopename[4]="sc46";
    isotopename[5]="v48";
*/
    Char_t *samplename = "PCBXX";

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    gROOT->LoadMacro("./get_th1f.C");

    for (Int_t iisotope = 0; iisotope < nisotopes; iisotope++) {

        cout << " now produce " << isotopename[iisotope] << " smeared MC " << endl;

//-----inputfilename need modification
        Char_t *inputfilename = Form("../%s/%s/rootfile/hMC_%s_%s_simplified.root", samplename, isotopename[iisotope],
                                     samplename, isotopename[iisotope]);

//-----outputfilename need modification
        Char_t *outputfilename = Form("../%s/%s/rootfile/hMC_%s_%s_smear.root", samplename, isotopename[iisotope],
                                      samplename, isotopename[iisotope]);

//------------------------------------------------------
// codes below needn't modifications
//------------------------------------------------------

        Int_t number_of_generated_events = 1.0; // unit million
        cout << " input file " << inputfilename << endl;
        cout << " output file " << outputfilename << endl;

//---> input ntuples
        TFile *finput = new TFile(inputfilename);
        TTree *t1 = (TTree *) gDirectory->Get("t1");
//Declaration of leaves types
        Float_t Cha_MCAEnergy;
        Float_t Cha_Energy;
        // Set branch addresses.
        t1->SetBranchAddress("Cha_MCAEnergy", &Cha_MCAEnergy);
        t1->SetBranchAddress("Cha_Energy", &Cha_Energy);
        Int_t nentries = t1->GetEntries();
        cout << " total number of events " << nentries << endl;
//---> output ntuples and histograms
        const Int_t MAXENTRY = 65536;
        Float_t escale = 0.0667788; // fake energy scale

        TH1F *hmca;
        hmca = new TH1F("hmca", "hmca",
                        MAXENTRY, 0.0, float(MAXENTRY));
        hmca->SetXTitle("ADC_counts");
        hmca->SetYTitle("Entries");

        TH1F *he;
        he = new TH1F("he", "he", 6000, 0.0, 3000.0);
        he->SetXTitle("E [keV]");
        he->SetYTitle("Entries / 0.5keV");

        TH1F *herun;
        herun = new TH1F("herun", "herun", 6000, 0.0, 3000.0);
        herun->SetXTitle("energy");
        herun->SetYTitle("Entries / keV day");

        TH1F *helive;
        helive = new TH1F("helive", "helive", 6000, 0.0, 3000.0);
        helive->SetXTitle("energy");
        helive->SetYTitle("Entries / keV day");

        TH1F *hparameter;
        hparameter = new TH1F("hparameter", "hparameter", 12, 0.0, 12.0);

        TFile *foutput = new TFile(outputfilename, "recreate");
        TTree *toutput = new TTree("t1", "Energy results");
        toutput->SetMaxTreeSize(100e9);
        toutput->Branch("Cha_MCAEnergy", &Cha_MCAEnergy, "Cha_MCAEnergy/F");
        toutput->Branch("Cha_Energy", &Cha_Energy, "Cha_Energy/F");

//---> start looping over ntuples
        Float_t sigma_e;
        TRandom3 r;
        for (Int_t ievt = 0; ievt < nentries; ievt++) {
            t1->GetEntry(ievt);
            if (ievt % 10000 == 0) cout << " now event " << ievt << endl;

            sigma_e = 0.512 + 0.0002632 * Cha_Energy;
            Cha_Energy += r.Gaus(0.0, sigma_e);
            foutput->cd();
            toutput->Fill();

            hmca->Fill(Cha_MCAEnergy); // MCAEnergy will not be smeared
            he->Fill(Cha_Energy);
            herun->Fill(Cha_Energy);
            helive->Fill(Cha_Energy);

        }
        toutput->Write();

        Float_t runtime = number_of_generated_events;
        Float_t livetime = runtime;
        hparameter->SetBinContent(1, runtime);
        hparameter->SetBinContent(2, livetime);
        hparameter->SetBinContent(12, escale);

        foutput->cd();
        hmca->Write();
        he->Write();
        herun->Write();
        helive->Write();
        hparameter->Write();

        foutput->Close();
        finput->Close();


        cout << " now finish " << isotopename[iisotope] << " smeared MC " << endl;
        cout << "=============================================================" << endl;
    }
}


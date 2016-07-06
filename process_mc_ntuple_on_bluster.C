/*
  this macro converts the MC ntuple
  to the same format as the data file 
  namely, ntuple varialbes Cha_MCAEnergy, Cha_Energy.

  histograms:
  hmca
  he
  herun
  helive

  hparameter
  1: run time
  2: live time
  3: rise time
  4: flat top
  5: decay time
  6: dynamic range
  7: trigger threshold
  8: input count rate
  9: output count rate
  10: FTDT
  11: pile up content
  12: energy scale
*/


void process_mc_ntuple_on_bluster()
{

//--------------------------------------------------
// following variables need user input
//--------------------------------------------------
    const Int_t nisotopes = 7;
    Char_t *isotopename[nisotopes];

    isotopename[0] = "co60";
    isotopename[1] = "cs137";
    isotopename[2] = "k40";
    isotopename[3] = "u238";
    isotopename[4] = "u235";
    isotopename[5] = "th232";
    isotopename[6] = "pb210";

    Char_t *samplename = "PCBXX";

//Reset ROOT and connect tree file
    gROOT->Reset();
    gROOT->SetStyle("Plain");

    for (Int_t iisotope = 0; iisotope < nisotopes; iisotope++) {
        cout << " now produce " << isotopename[iisotope] << " MC " << endl;

//-----inputfilename need modification
        Char_t *inputfilename = Form("../%s/%s/rootfile/hMC_%s_%s.root", samplename, isotopename[iisotope], samplename,
                                     isotopename[iisotope]);

//-----outputfilename need modification
        Char_t *outputfilename = Form("../%s/%s/rootfile/hMC_%s_%s_simplified.root", samplename, isotopename[iisotope],
                                      samplename, isotopename[iisotope]);

        cout << " input file " << inputfilename << endl;
        cout << " output file " << outputfilename << endl;

//-----this parameter need to be modified in some cases
        Float_t number_of_generated_events = 1.0; // unit Million
        Float_t weight = 1.0; // unit kg (or piece etc)
        Float_t contamination_level = 1.0;  // Bq/kg

//-------------------------------------------
// output histograms and ntuple
//-------------------------------------------
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

        Float_t Cha_MCAEnergy;
        Float_t Cha_Energy;
//     cout<<"here defined"<<endl;
        TFile *foutput = new TFile(outputfilename, "recreate");
        TTree *toutput = new TTree("t1", "Energy results");
//    toutput->SetMaxTreeSize(100e9);
        toutput->Branch("Cha_MCAEnergy", &Cha_MCAEnergy, "Cha_MCAEnergy/F");
        toutput->Branch("Cha_Energy", &Cha_Energy, "Cha_Energy/F");

//-------------------------------------------
// input ntuple
//-------------------------------------------

        TChain *t1 = new TChain("t1");
        t1->AddFile(inputfilename);
//      cout<<"here again"<<endl;
//Declaration of leaves types
        Int_t eventnumber;
        Int_t vertex_totnum;

        Int_t hits_totnum;
        Float_t hits_tote;

        t1->SetBranchAddress("hits_totnum", &hits_totnum);
        t1->SetBranchAddress("hits_tote", &hits_tote);


        Int_t nentries = t1->GetEntries();
        cout << " total number of MC events " << nentries << endl;

        for (Int_t ievt = 0; ievt < nentries; ievt++) {
            t1->GetEntry(ievt);
            if (ievt % 10000 == 0) cout << " now event " << ievt << endl;

            Cha_Energy = hits_tote;
            Cha_MCAEnergy = hits_tote / escale;
            foutput->cd();
            toutput->Fill();

            hmca->Fill(Cha_MCAEnergy);
            he->Fill(Cha_Energy);
            herun->Fill(Cha_Energy);
            helive->Fill(Cha_Energy);
        }
        toutput->Write();

        Float_t number_of_decays_per_day = weight * contamination_level * 24.0 * 3600.0 / 1000000.0; // unit M
        Float_t scale_factor = number_of_decays_per_day / 2.0 / number_of_generated_events;

        herun->Scale(scale_factor);  // herun y axis unit events/day/keV
        helive->Scale(scale_factor);

        Float_t runtime = number_of_generated_events / number_of_decays_per_day * 24.0 * 3600.0; //unit second
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

        cout << " now finish " << isotopename[iisotope] << " MC " << endl;
        cout << "=============================================================" << endl;
    }
}

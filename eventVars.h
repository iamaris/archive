//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Aug 25 16:24:09 2013 by ROOT version 5.99/02
// from TTree eventVars/Event variables
// found on file: A2012_minintuple_1.root
//////////////////////////////////////////////////////////

#ifndef eventVars_h
#define eventVars_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class eventVars {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         met;
   Float_t         metPhi;
   Float_t         mht;
   Float_t         mhtPhi;
   Float_t         mtElectron;
   Float_t         mtElectronPhoton;
   Float_t         mtMuon;
   Float_t         mtMuonPhoton;
   Float_t         enuMomentum;
   Float_t         munuMomentum;
   Bool_t          HLT_Ele27_WP80;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoMu24_eta2p1;
   Bool_t          HLT_Mu22_Photon22_CaloIdL;
   Bool_t          HLT_Photon135;
   Bool_t          HLT_Photon150;
   Bool_t          HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60;
   Bool_t          HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60;
   Bool_t          HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50;
   Bool_t          HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50;
   Bool_t          HLT_Photon60_CaloIdL_HT300;
   Bool_t          HLT_Photon60_CaloIdL_MHT70;
   Bool_t          HLT_Photon70_CaloIdXL_PFHT400;
   Bool_t          HLT_Photon70_CaloIdXL_PFHT500;
   Bool_t          HLT_Photon70_CaloIdXL_PFMET100;
   Bool_t          HLT_Photon70_CaloIdXL_PFNoPUHT400;
   Bool_t          HLT_Photon70_CaloIdXL_PFNoPUHT500;
   Float_t         pthat;
   Float_t         genMet;
   Float_t         genMetPhi;
   UInt_t          gen_size;
   UShort_t        gen_status[1];   //[gen.size]
   Short_t         gen_charge[1];   //[gen.size]
   Short_t         gen_motherIndex[1];   //[gen.size]
   Int_t           gen_pdgId[1];   //[gen.size]
   Float_t         gen_vx[1];   //[gen.size]
   Float_t         gen_vy[1];   //[gen.size]
   Float_t         gen_vz[1];   //[gen.size]
   Float_t         gen_pt[1];   //[gen.size]
   Float_t         gen_eta[1];   //[gen.size]
   Float_t         gen_phi[1];   //[gen.size]
   Float_t         gen_mass[1];   //[gen.size]
   Float_t         gen_px[1];   //[gen.size]
   Float_t         gen_py[1];   //[gen.size]
   Float_t         gen_pz[1];   //[gen.size]
   Float_t         gen_energy[1];   //[gen.size]

   // List of branches
   TBranch        *b_met;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_mht;   //!
   TBranch        *b_mhtPhi;   //!
   TBranch        *b_mtElectron;   //!
   TBranch        *b_mtElectronPhoton;   //!
   TBranch        *b_mtMuon;   //!
   TBranch        *b_mtMuonPhoton;   //!
   TBranch        *b_enuMomentum;   //!
   TBranch        *b_munuMomentum;   //!
   TBranch        *b_HLT_Ele27_WP80;   //!
   TBranch        *b_HLT_IsoMu24;   //!
   TBranch        *b_HLT_IsoMu24_eta2p1;   //!
   TBranch        *b_HLT_Mu22_Photon22_CaloIdL;   //!
   TBranch        *b_HLT_Photon135;   //!
   TBranch        *b_HLT_Photon150;   //!
   TBranch        *b_HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60;   //!
   TBranch        *b_HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60;   //!
   TBranch        *b_HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50;   //!
   TBranch        *b_HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50;   //!
   TBranch        *b_HLT_Photon60_CaloIdL_HT300;   //!
   TBranch        *b_HLT_Photon60_CaloIdL_MHT70;   //!
   TBranch        *b_HLT_Photon70_CaloIdXL_PFHT400;   //!
   TBranch        *b_HLT_Photon70_CaloIdXL_PFHT500;   //!
   TBranch        *b_HLT_Photon70_CaloIdXL_PFMET100;   //!
   TBranch        *b_HLT_Photon70_CaloIdXL_PFNoPUHT400;   //!
   TBranch        *b_HLT_Photon70_CaloIdXL_PFNoPUHT500;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_genMet;   //!
   TBranch        *b_genMetPhi;   //!
   TBranch        *b_gen_size;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_charge;   //!
   TBranch        *b_gen_motherIndex;   //!
   TBranch        *b_gen_pdgId;   //!
   TBranch        *b_gen_vx;   //!
   TBranch        *b_gen_vy;   //!
   TBranch        *b_gen_vz;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_mass;   //!
   TBranch        *b_gen_px;   //!
   TBranch        *b_gen_py;   //!
   TBranch        *b_gen_pz;   //!
   TBranch        *b_gen_energy;   //!

   eventVars(TTree *tree=0);
   virtual ~eventVars();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef eventVars_cxx
eventVars::eventVars(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("A2012_minintuple_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("A2012_minintuple_1.root");
      }
      f->GetObject("eventVars",tree);

   }
   Init(tree);
}

eventVars::~eventVars()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t eventVars::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t eventVars::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void eventVars::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("mht", &mht, &b_mht);
   fChain->SetBranchAddress("mhtPhi", &mhtPhi, &b_mhtPhi);
   fChain->SetBranchAddress("mtElectron", &mtElectron, &b_mtElectron);
   fChain->SetBranchAddress("mtElectronPhoton", &mtElectronPhoton, &b_mtElectronPhoton);
   fChain->SetBranchAddress("mtMuon", &mtMuon, &b_mtMuon);
   fChain->SetBranchAddress("mtMuonPhoton", &mtMuonPhoton, &b_mtMuonPhoton);
   fChain->SetBranchAddress("enuMomentum", &enuMomentum, &b_enuMomentum);
   fChain->SetBranchAddress("munuMomentum", &munuMomentum, &b_munuMomentum);
   fChain->SetBranchAddress("HLT_Ele27_WP80", &HLT_Ele27_WP80, &b_HLT_Ele27_WP80);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu24_eta2p1", &HLT_IsoMu24_eta2p1, &b_HLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("HLT_Mu22_Photon22_CaloIdL", &HLT_Mu22_Photon22_CaloIdL, &b_HLT_Mu22_Photon22_CaloIdL);
   fChain->SetBranchAddress("HLT_Photon135", &HLT_Photon135, &b_HLT_Photon135);
   fChain->SetBranchAddress("HLT_Photon150", &HLT_Photon150, &b_HLT_Photon150);
   fChain->SetBranchAddress("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60", &HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60, &b_HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60);
   fChain->SetBranchAddress("HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60", &HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60, &b_HLT_Photon26_R9Id85_IR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60);
   fChain->SetBranchAddress("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50", &HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50, &b_HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50);
   fChain->SetBranchAddress("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50", &HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50, &b_HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50);
   fChain->SetBranchAddress("HLT_Photon60_CaloIdL_HT300", &HLT_Photon60_CaloIdL_HT300, &b_HLT_Photon60_CaloIdL_HT300);
   fChain->SetBranchAddress("HLT_Photon60_CaloIdL_MHT70", &HLT_Photon60_CaloIdL_MHT70, &b_HLT_Photon60_CaloIdL_MHT70);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdXL_PFHT400", &HLT_Photon70_CaloIdXL_PFHT400, &b_HLT_Photon70_CaloIdXL_PFHT400);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdXL_PFHT500", &HLT_Photon70_CaloIdXL_PFHT500, &b_HLT_Photon70_CaloIdXL_PFHT500);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdXL_PFMET100", &HLT_Photon70_CaloIdXL_PFMET100, &b_HLT_Photon70_CaloIdXL_PFMET100);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdXL_PFNoPUHT400", &HLT_Photon70_CaloIdXL_PFNoPUHT400, &b_HLT_Photon70_CaloIdXL_PFNoPUHT400);
   fChain->SetBranchAddress("HLT_Photon70_CaloIdXL_PFNoPUHT500", &HLT_Photon70_CaloIdXL_PFNoPUHT500, &b_HLT_Photon70_CaloIdXL_PFNoPUHT500);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("genMet", &genMet, &b_genMet);
   fChain->SetBranchAddress("genMetPhi", &genMetPhi, &b_genMetPhi);
   fChain->SetBranchAddress("gen.size", &gen_size, &b_gen_size);
   fChain->SetBranchAddress("gen.status", &gen_status, &b_gen_status);
   fChain->SetBranchAddress("gen.charge", &gen_charge, &b_gen_charge);
   fChain->SetBranchAddress("gen.motherIndex", &gen_motherIndex, &b_gen_motherIndex);
   fChain->SetBranchAddress("gen.pdgId", &gen_pdgId, &b_gen_pdgId);
   fChain->SetBranchAddress("gen.vx", &gen_vx, &b_gen_vx);
   fChain->SetBranchAddress("gen.vy", &gen_vy, &b_gen_vy);
   fChain->SetBranchAddress("gen.vz", &gen_vz, &b_gen_vz);
   fChain->SetBranchAddress("gen.pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen.eta", &gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen.phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen.mass", &gen_mass, &b_gen_mass);
   fChain->SetBranchAddress("gen.px", &gen_px, &b_gen_px);
   fChain->SetBranchAddress("gen.py", &gen_py, &b_gen_py);
   fChain->SetBranchAddress("gen.pz", &gen_pz, &b_gen_pz);
   fChain->SetBranchAddress("gen.energy", &gen_energy, &b_gen_energy);
   Notify();
}

Bool_t eventVars::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void eventVars::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t eventVars::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef eventVars_cxx

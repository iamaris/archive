//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Aug 25 16:24:46 2013 by ROOT version 5.99/02
// from TTree allObjects/All objects
// found on file: A2012_minintuple_1.root
//////////////////////////////////////////////////////////

#ifndef allObjects_h
#define allObjects_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class allObjects {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          runNumber;
   UInt_t          lumiNumber;
   UInt_t          eventNumber;
   UInt_t          photon_size;
   Float_t         photon_pt[11];   //[photon.size]
   Float_t         photon_eta[11];   //[photon.size]
   Float_t         photon_phi[11];   //[photon.size]
   Float_t         photon_px[11];   //[photon.size]
   Float_t         photon_py[11];   //[photon.size]
   Float_t         photon_pz[11];   //[photon.size]
   Float_t         photon_energy[11];   //[photon.size]
   Float_t         photon_hOverE[11];   //[photon.size]
   Float_t         photon_sigmaIetaIeta[11];   //[photon.size]
   Float_t         photon_sigmaIphiIphi[11];   //[photon.size]
   Float_t         photon_etaWidth[11];   //[photon.size]
   Float_t         photon_phiWidth[11];   //[photon.size]
   Float_t         photon_r9[11];   //[photon.size]
   Float_t         photon_r5[11];   //[photon.size]
   Float_t         photon_trackerIso[11];   //[photon.size]
   Float_t         photon_ecalIso[11];   //[photon.size]
   Float_t         photon_hcalIso[11];   //[photon.size]
   Float_t         photon_chargedHadronIso[11];   //[photon.size]
   Float_t         photon_neutralHadronIso[11];   //[photon.size]
   Float_t         photon_photonIso[11];   //[photon.size]
   Short_t         photon_iSubdet[11];   //[photon.size]
   UChar_t         photon_nPixelSeeds[11];   //[photon.size]
   UChar_t         photon_nClusters[11];   //[photon.size]
   Bool_t          photon_hasMatchedElectron[11];   //[photon.size]
   Bool_t          photon_electronVetoBit[11];   //[photon.size]
   Bool_t          photon_looseElectronVetoBit[11];   //[photon.size]
   Bool_t          photon_isLoose[11];   //[photon.size]
   Bool_t          photon_isMedium[11];   //[photon.size]
   Bool_t          photon_isTight[11];   //[photon.size]
   Bool_t          photon_isLoosePix[11];   //[photon.size]
   Bool_t          photon_isMediumPix[11];   //[photon.size]
   Bool_t          photon_isTightPix[11];   //[photon.size]
   Bool_t          photon_isLooseLV[11];   //[photon.size]
   Bool_t          photon_isMediumLV[11];   //[photon.size]
   Bool_t          photon_isTightLV[11];   //[photon.size]
   UInt_t          electron_size;
   Float_t         electron_pt[8];   //[electron.size]
   Float_t         electron_eta[8];   //[electron.size]
   Float_t         electron_phi[8];   //[electron.size]
   Float_t         electron_px[8];   //[electron.size]
   Float_t         electron_py[8];   //[electron.size]
   Float_t         electron_pz[8];   //[electron.size]
   Float_t         electron_energy[8];   //[electron.size]
   Float_t         electron_combRelSubdetIso[8];   //[electron.size]
   Float_t         electron_combRelIso[8];   //[electron.size]
   Float_t         electron_deltaEta[8];   //[electron.size]
   Float_t         electron_deltaPhi[8];   //[electron.size]
   Float_t         electron_sigmaIetaIeta[8];   //[electron.size]
   Float_t         electron_sigmaIphiIphi[8];   //[electron.size]
   Float_t         electron_r9[8];   //[electron.size]
   Float_t         electron_r5[8];   //[electron.size]
   Float_t         electron_etaWidth[8];   //[electron.size]
   Float_t         electron_phiWidth[8];   //[electron.size]
   Float_t         electron_hOverE[8];   //[electron.size]
   Float_t         electron_d0[8];   //[electron.size]
   Float_t         electron_dz[8];   //[electron.size]
   Float_t         electron_epDiff[8];   //[electron.size]
   Float_t         electron_vtxFitProb[8];   //[electron.size]
   Float_t         electron_dCot[8];   //[electron.size]
   Float_t         electron_dist[8];   //[electron.size]
   Short_t         electron_iSubdet[8];   //[electron.size]
   UChar_t         electron_nClusters[8];   //[electron.size]
   UChar_t         electron_nPixelHits[8];   //[electron.size]
   UChar_t         electron_nMissingHits[8];   //[electron.size]
   Bool_t          electron_passConversionVeto[8];   //[electron.size]
   Bool_t          electron_isVeto[8];   //[electron.size]
   Bool_t          electron_isLoose[8];   //[electron.size]
   Bool_t          electron_isMedium[8];   //[electron.size]
   Bool_t          electron_isTight[8];   //[electron.size]
   UInt_t          muon_size;
   Float_t         muon_pt[27];   //[muon.size]
   Float_t         muon_eta[27];   //[muon.size]
   Float_t         muon_phi[27];   //[muon.size]
   Float_t         muon_px[27];   //[muon.size]
   Float_t         muon_py[27];   //[muon.size]
   Float_t         muon_pz[27];   //[muon.size]
   Float_t         muon_energy[27];   //[muon.size]
   Float_t         muon_normChi2[27];   //[muon.size]
   Float_t         muon_dxy[27];   //[muon.size]
   Float_t         muon_dz[27];   //[muon.size]
   Float_t         muon_combRelSubdetIso[27];   //[muon.size]
   Float_t         muon_combRelIso[27];   //[muon.size]
   Short_t         muon_iSubdet[27];   //[muon.size]
   UChar_t         muon_nMatchedStations[27];   //[muon.size]
   UChar_t         muon_nLayersWithMmt[27];   //[muon.size]
   UChar_t         muon_nValidMuonHits[27];   //[muon.size]
   UChar_t         muon_nValidPixelHits[27];   //[muon.size]
   Bool_t          muon_isGlobalMuon[27];   //[muon.size]
   Bool_t          muon_isPFMuon[27];   //[muon.size]
   Bool_t          muon_hasInnerTrack[27];   //[muon.size]
   Bool_t          muon_hasGlobalTrack[27];   //[muon.size]
   Bool_t          muon_hasBestTrack[27];   //[muon.size]
   Bool_t          muon_isLoose[27];   //[muon.size]
   Bool_t          muon_isTight[27];   //[muon.size]
   UInt_t          jet_size;
   Float_t         jet_pt[18];   //[jet.size]
   Float_t         jet_eta[18];   //[jet.size]
   Float_t         jet_phi[18];   //[jet.size]
   Float_t         jet_px[18];   //[jet.size]
   Float_t         jet_py[18];   //[jet.size]
   Float_t         jet_pz[18];   //[jet.size]
   Float_t         jet_energy[18];   //[jet.size]
   Float_t         jet_jecScale[18];   //[jet.size]
   Float_t         jet_chFraction[18];   //[jet.size]
   Float_t         jet_nhFraction[18];   //[jet.size]
   Float_t         jet_ceFraction[18];   //[jet.size]
   Float_t         jet_neFraction[18];   //[jet.size]
   Short_t         jet_iSubdet[18];   //[jet.size]
   UChar_t         jet_nConstituents[18];   //[jet.size]
   UChar_t         jet_nCharged[18];   //[jet.size]
   Bool_t          jet_isLoose[18];   //[jet.size]
   UInt_t          vertex_size;
   Float_t         vertex_x[38];   //[vertex.size]
   Float_t         vertex_y[38];   //[vertex.size]
   Float_t         vertex_z[38];   //[vertex.size]
   Float_t         vertex_rho[38];   //[vertex.size]
   Float_t         vertex_sumPt2[38];   //[vertex.size]
   Float_t         vertex_chi2[38];   //[vertex.size]
   Float_t         vertex_ndof[38];   //[vertex.size]
   Bool_t          vertex_isGood[38];   //[vertex.size]
   Float_t         photon_dRGen[11];   //[photon.size]
   Float_t         photon_genIso[11];   //[photon.size]
   Int_t           photon_nearestGen[11];   //[photon.size]
   Float_t         photon_dRJet[11];   //[photon.size]
   Float_t         photon_dRNextJet[11];   //[photon.size]
   Float_t         photon_dRPF[11];   //[photon.size]
   Short_t         photon_nearestPF[11];   //[photon.size]
   Bool_t          photon_pfIsPU[11];   //[photon.size]
   Float_t         electron_dRGen[8];   //[electron.size]
   Float_t         electron_genIso[8];   //[electron.size]
   Int_t           electron_nearestGen[8];   //[electron.size]
   Float_t         electron_dRJet[8];   //[electron.size]
   Float_t         electron_dRNextJet[8];   //[electron.size]
   Float_t         electron_dRPhoton[8];   //[electron.size]
   Float_t         electron_dRNextPhoton[8];   //[electron.size]
   Float_t         electron_dRPF[8];   //[electron.size]
   Short_t         electron_nearestPF[8];   //[electron.size]
   Bool_t          electron_pfIsPU[8];   //[electron.size]
   Float_t         muon_dRGen[27];   //[muon.size]
   Float_t         muon_genIso[27];   //[muon.size]
   Int_t           muon_nearestGen[27];   //[muon.size]
   Float_t         muon_dRJet[27];   //[muon.size]
   Float_t         muon_dRNextJet[27];   //[muon.size]
   Float_t         muon_dRPhoton[27];   //[muon.size]
   Float_t         muon_dRNextPhoton[27];   //[muon.size]
   Float_t         muon_dRPF[27];   //[muon.size]
   Short_t         muon_nearestPF[27];   //[muon.size]
   Bool_t          muon_pfIsPU[27];   //[muon.size]
   Float_t         jet_dRGen[18];   //[jet.size]
   Float_t         jet_genSumPt[18];   //[jet.size]
   Int_t           jet_nearestGen[18];   //[jet.size]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_lumiNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_photon_size;   //!
   TBranch        *b_photon_pt;   //!
   TBranch        *b_photon_eta;   //!
   TBranch        *b_photon_phi;   //!
   TBranch        *b_photon_px;   //!
   TBranch        *b_photon_py;   //!
   TBranch        *b_photon_pz;   //!
   TBranch        *b_photon_energy;   //!
   TBranch        *b_photon_hOverE;   //!
   TBranch        *b_photon_sigmaIetaIeta;   //!
   TBranch        *b_photon_sigmaIphiIphi;   //!
   TBranch        *b_photon_etaWidth;   //!
   TBranch        *b_photon_phiWidth;   //!
   TBranch        *b_photon_r9;   //!
   TBranch        *b_photon_r5;   //!
   TBranch        *b_photon_trackerIso;   //!
   TBranch        *b_photon_ecalIso;   //!
   TBranch        *b_photon_hcalIso;   //!
   TBranch        *b_photon_chargedHadronIso;   //!
   TBranch        *b_photon_neutralHadronIso;   //!
   TBranch        *b_photon_photonIso;   //!
   TBranch        *b_photon_iSubdet;   //!
   TBranch        *b_photon_nPixelSeeds;   //!
   TBranch        *b_photon_nClusters;   //!
   TBranch        *b_photon_hasMatchedElectron;   //!
   TBranch        *b_photon_electronVetoBit;   //!
   TBranch        *b_photon_looseElectronVetoBit;   //!
   TBranch        *b_photon_isLoose;   //!
   TBranch        *b_photon_isMedium;   //!
   TBranch        *b_photon_isTight;   //!
   TBranch        *b_photon_isLoosePix;   //!
   TBranch        *b_photon_isMediumPix;   //!
   TBranch        *b_photon_isTightPix;   //!
   TBranch        *b_photon_isLooseLV;   //!
   TBranch        *b_photon_isMediumLV;   //!
   TBranch        *b_photon_isTightLV;   //!
   TBranch        *b_electron_size;   //!
   TBranch        *b_electron_pt;   //!
   TBranch        *b_electron_eta;   //!
   TBranch        *b_electron_phi;   //!
   TBranch        *b_electron_px;   //!
   TBranch        *b_electron_py;   //!
   TBranch        *b_electron_pz;   //!
   TBranch        *b_electron_energy;   //!
   TBranch        *b_electron_combRelSubdetIso;   //!
   TBranch        *b_electron_combRelIso;   //!
   TBranch        *b_electron_deltaEta;   //!
   TBranch        *b_electron_deltaPhi;   //!
   TBranch        *b_electron_sigmaIetaIeta;   //!
   TBranch        *b_electron_sigmaIphiIphi;   //!
   TBranch        *b_electron_r9;   //!
   TBranch        *b_electron_r5;   //!
   TBranch        *b_electron_etaWidth;   //!
   TBranch        *b_electron_phiWidth;   //!
   TBranch        *b_electron_hOverE;   //!
   TBranch        *b_electron_d0;   //!
   TBranch        *b_electron_dz;   //!
   TBranch        *b_electron_epDiff;   //!
   TBranch        *b_electron_vtxFitProb;   //!
   TBranch        *b_electron_dCot;   //!
   TBranch        *b_electron_dist;   //!
   TBranch        *b_electron_iSubdet;   //!
   TBranch        *b_electron_nClusters;   //!
   TBranch        *b_electron_nPixelHits;   //!
   TBranch        *b_electron_nMissingHits;   //!
   TBranch        *b_electron_passConversionVeto;   //!
   TBranch        *b_electron_isVeto;   //!
   TBranch        *b_electron_isLoose;   //!
   TBranch        *b_electron_isMedium;   //!
   TBranch        *b_electron_isTight;   //!
   TBranch        *b_muon_size;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_px;   //!
   TBranch        *b_muon_py;   //!
   TBranch        *b_muon_pz;   //!
   TBranch        *b_muon_energy;   //!
   TBranch        *b_muon_normChi2;   //!
   TBranch        *b_muon_dxy;   //!
   TBranch        *b_muon_dz;   //!
   TBranch        *b_muon_combRelSubdetIso;   //!
   TBranch        *b_muon_combRelIso;   //!
   TBranch        *b_muon_iSubdet;   //!
   TBranch        *b_muon_nMatchedStations;   //!
   TBranch        *b_muon_nLayersWithMmt;   //!
   TBranch        *b_muon_nValidMuonHits;   //!
   TBranch        *b_muon_nValidPixelHits;   //!
   TBranch        *b_muon_isGlobalMuon;   //!
   TBranch        *b_muon_isPFMuon;   //!
   TBranch        *b_muon_hasInnerTrack;   //!
   TBranch        *b_muon_hasGlobalTrack;   //!
   TBranch        *b_muon_hasBestTrack;   //!
   TBranch        *b_muon_isLoose;   //!
   TBranch        *b_muon_isTight;   //!
   TBranch        *b_jet_size;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_px;   //!
   TBranch        *b_jet_py;   //!
   TBranch        *b_jet_pz;   //!
   TBranch        *b_jet_energy;   //!
   TBranch        *b_jet_jecScale;   //!
   TBranch        *b_jet_chFraction;   //!
   TBranch        *b_jet_nhFraction;   //!
   TBranch        *b_jet_ceFraction;   //!
   TBranch        *b_jet_neFraction;   //!
   TBranch        *b_jet_iSubdet;   //!
   TBranch        *b_jet_nConstituents;   //!
   TBranch        *b_jet_nCharged;   //!
   TBranch        *b_jet_isLoose;   //!
   TBranch        *b_vertex_size;   //!
   TBranch        *b_vertex_x;   //!
   TBranch        *b_vertex_y;   //!
   TBranch        *b_vertex_z;   //!
   TBranch        *b_vertex_rho;   //!
   TBranch        *b_vertex_sumPt2;   //!
   TBranch        *b_vertex_chi2;   //!
   TBranch        *b_vertex_ndof;   //!
   TBranch        *b_vertex_isGood;   //!
   TBranch        *b_photon_dRGen;   //!
   TBranch        *b_photon_genIso;   //!
   TBranch        *b_photon_nearestGen;   //!
   TBranch        *b_photon_dRJet;   //!
   TBranch        *b_photon_dRNextJet;   //!
   TBranch        *b_photon_dRPF;   //!
   TBranch        *b_photon_nearestPF;   //!
   TBranch        *b_photon_pfIsPU;   //!
   TBranch        *b_electron_dRGen;   //!
   TBranch        *b_electron_genIso;   //!
   TBranch        *b_electron_nearestGen;   //!
   TBranch        *b_electron_dRJet;   //!
   TBranch        *b_electron_dRNextJet;   //!
   TBranch        *b_electron_dRPhoton;   //!
   TBranch        *b_electron_dRNextPhoton;   //!
   TBranch        *b_electron_dRPF;   //!
   TBranch        *b_electron_nearestPF;   //!
   TBranch        *b_electron_pfIsPU;   //!
   TBranch        *b_muon_dRGen;   //!
   TBranch        *b_muon_genIso;   //!
   TBranch        *b_muon_nearestGen;   //!
   TBranch        *b_muon_dRJet;   //!
   TBranch        *b_muon_dRNextJet;   //!
   TBranch        *b_muon_dRPhoton;   //!
   TBranch        *b_muon_dRNextPhoton;   //!
   TBranch        *b_muon_dRPF;   //!
   TBranch        *b_muon_nearestPF;   //!
   TBranch        *b_muon_pfIsPU;   //!
   TBranch        *b_jet_dRGen;   //!
   TBranch        *b_jet_genSumPt;   //!
   TBranch        *b_jet_nearestGen;   //!

   allObjects(TTree *tree=0);
   virtual ~allObjects();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef allObjects_cxx
allObjects::allObjects(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("A2012_minintuple_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("A2012_minintuple_1.root");
      }
      f->GetObject("allObjects",tree);

   }
   Init(tree);
}

allObjects::~allObjects()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t allObjects::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t allObjects::LoadTree(Long64_t entry)
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

void allObjects::Init(TTree *tree)
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

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("lumiNumber", &lumiNumber, &b_lumiNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("photon.size", &photon_size, &b_photon_size);
   fChain->SetBranchAddress("photon.pt", photon_pt, &b_photon_pt);
   fChain->SetBranchAddress("photon.eta", photon_eta, &b_photon_eta);
   fChain->SetBranchAddress("photon.phi", photon_phi, &b_photon_phi);
   fChain->SetBranchAddress("photon.px", photon_px, &b_photon_px);
   fChain->SetBranchAddress("photon.py", photon_py, &b_photon_py);
   fChain->SetBranchAddress("photon.pz", photon_pz, &b_photon_pz);
   fChain->SetBranchAddress("photon.energy", photon_energy, &b_photon_energy);
   fChain->SetBranchAddress("photon.hOverE", photon_hOverE, &b_photon_hOverE);
   fChain->SetBranchAddress("photon.sigmaIetaIeta", photon_sigmaIetaIeta, &b_photon_sigmaIetaIeta);
   fChain->SetBranchAddress("photon.sigmaIphiIphi", photon_sigmaIphiIphi, &b_photon_sigmaIphiIphi);
   fChain->SetBranchAddress("photon.etaWidth", photon_etaWidth, &b_photon_etaWidth);
   fChain->SetBranchAddress("photon.phiWidth", photon_phiWidth, &b_photon_phiWidth);
   fChain->SetBranchAddress("photon.r9", photon_r9, &b_photon_r9);
   fChain->SetBranchAddress("photon.r5", photon_r5, &b_photon_r5);
   fChain->SetBranchAddress("photon.trackerIso", photon_trackerIso, &b_photon_trackerIso);
   fChain->SetBranchAddress("photon.ecalIso", photon_ecalIso, &b_photon_ecalIso);
   fChain->SetBranchAddress("photon.hcalIso", photon_hcalIso, &b_photon_hcalIso);
   fChain->SetBranchAddress("photon.chargedHadronIso", photon_chargedHadronIso, &b_photon_chargedHadronIso);
   fChain->SetBranchAddress("photon.neutralHadronIso", photon_neutralHadronIso, &b_photon_neutralHadronIso);
   fChain->SetBranchAddress("photon.photonIso", photon_photonIso, &b_photon_photonIso);
   fChain->SetBranchAddress("photon.iSubdet", photon_iSubdet, &b_photon_iSubdet);
   fChain->SetBranchAddress("photon.nPixelSeeds", photon_nPixelSeeds, &b_photon_nPixelSeeds);
   fChain->SetBranchAddress("photon.nClusters", photon_nClusters, &b_photon_nClusters);
   fChain->SetBranchAddress("photon.hasMatchedElectron", photon_hasMatchedElectron, &b_photon_hasMatchedElectron);
   fChain->SetBranchAddress("photon.electronVetoBit", photon_electronVetoBit, &b_photon_electronVetoBit);
   fChain->SetBranchAddress("photon.looseElectronVetoBit", photon_looseElectronVetoBit, &b_photon_looseElectronVetoBit);
   fChain->SetBranchAddress("photon.isLoose", photon_isLoose, &b_photon_isLoose);
   fChain->SetBranchAddress("photon.isMedium", photon_isMedium, &b_photon_isMedium);
   fChain->SetBranchAddress("photon.isTight", photon_isTight, &b_photon_isTight);
   fChain->SetBranchAddress("photon.isLoosePix", photon_isLoosePix, &b_photon_isLoosePix);
   fChain->SetBranchAddress("photon.isMediumPix", photon_isMediumPix, &b_photon_isMediumPix);
   fChain->SetBranchAddress("photon.isTightPix", photon_isTightPix, &b_photon_isTightPix);
   fChain->SetBranchAddress("photon.isLooseLV", photon_isLooseLV, &b_photon_isLooseLV);
   fChain->SetBranchAddress("photon.isMediumLV", photon_isMediumLV, &b_photon_isMediumLV);
   fChain->SetBranchAddress("photon.isTightLV", photon_isTightLV, &b_photon_isTightLV);
   fChain->SetBranchAddress("electron.size", &electron_size, &b_electron_size);
   fChain->SetBranchAddress("electron.pt", electron_pt, &b_electron_pt);
   fChain->SetBranchAddress("electron.eta", electron_eta, &b_electron_eta);
   fChain->SetBranchAddress("electron.phi", electron_phi, &b_electron_phi);
   fChain->SetBranchAddress("electron.px", electron_px, &b_electron_px);
   fChain->SetBranchAddress("electron.py", electron_py, &b_electron_py);
   fChain->SetBranchAddress("electron.pz", electron_pz, &b_electron_pz);
   fChain->SetBranchAddress("electron.energy", electron_energy, &b_electron_energy);
   fChain->SetBranchAddress("electron.combRelSubdetIso", electron_combRelSubdetIso, &b_electron_combRelSubdetIso);
   fChain->SetBranchAddress("electron.combRelIso", electron_combRelIso, &b_electron_combRelIso);
   fChain->SetBranchAddress("electron.deltaEta", electron_deltaEta, &b_electron_deltaEta);
   fChain->SetBranchAddress("electron.deltaPhi", electron_deltaPhi, &b_electron_deltaPhi);
   fChain->SetBranchAddress("electron.sigmaIetaIeta", electron_sigmaIetaIeta, &b_electron_sigmaIetaIeta);
   fChain->SetBranchAddress("electron.sigmaIphiIphi", electron_sigmaIphiIphi, &b_electron_sigmaIphiIphi);
   fChain->SetBranchAddress("electron.r9", electron_r9, &b_electron_r9);
   fChain->SetBranchAddress("electron.r5", electron_r5, &b_electron_r5);
   fChain->SetBranchAddress("electron.etaWidth", electron_etaWidth, &b_electron_etaWidth);
   fChain->SetBranchAddress("electron.phiWidth", electron_phiWidth, &b_electron_phiWidth);
   fChain->SetBranchAddress("electron.hOverE", electron_hOverE, &b_electron_hOverE);
   fChain->SetBranchAddress("electron.d0", electron_d0, &b_electron_d0);
   fChain->SetBranchAddress("electron.dz", electron_dz, &b_electron_dz);
   fChain->SetBranchAddress("electron.epDiff", electron_epDiff, &b_electron_epDiff);
   fChain->SetBranchAddress("electron.vtxFitProb", electron_vtxFitProb, &b_electron_vtxFitProb);
   fChain->SetBranchAddress("electron.dCot", electron_dCot, &b_electron_dCot);
   fChain->SetBranchAddress("electron.dist", electron_dist, &b_electron_dist);
   fChain->SetBranchAddress("electron.iSubdet", electron_iSubdet, &b_electron_iSubdet);
   fChain->SetBranchAddress("electron.nClusters", electron_nClusters, &b_electron_nClusters);
   fChain->SetBranchAddress("electron.nPixelHits", electron_nPixelHits, &b_electron_nPixelHits);
   fChain->SetBranchAddress("electron.nMissingHits", electron_nMissingHits, &b_electron_nMissingHits);
   fChain->SetBranchAddress("electron.passConversionVeto", electron_passConversionVeto, &b_electron_passConversionVeto);
   fChain->SetBranchAddress("electron.isVeto", electron_isVeto, &b_electron_isVeto);
   fChain->SetBranchAddress("electron.isLoose", electron_isLoose, &b_electron_isLoose);
   fChain->SetBranchAddress("electron.isMedium", electron_isMedium, &b_electron_isMedium);
   fChain->SetBranchAddress("electron.isTight", electron_isTight, &b_electron_isTight);
   fChain->SetBranchAddress("muon.size", &muon_size, &b_muon_size);
   fChain->SetBranchAddress("muon.pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon.eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon.phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon.px", muon_px, &b_muon_px);
   fChain->SetBranchAddress("muon.py", muon_py, &b_muon_py);
   fChain->SetBranchAddress("muon.pz", muon_pz, &b_muon_pz);
   fChain->SetBranchAddress("muon.energy", muon_energy, &b_muon_energy);
   fChain->SetBranchAddress("muon.normChi2", muon_normChi2, &b_muon_normChi2);
   fChain->SetBranchAddress("muon.dxy", muon_dxy, &b_muon_dxy);
   fChain->SetBranchAddress("muon.dz", muon_dz, &b_muon_dz);
   fChain->SetBranchAddress("muon.combRelSubdetIso", muon_combRelSubdetIso, &b_muon_combRelSubdetIso);
   fChain->SetBranchAddress("muon.combRelIso", muon_combRelIso, &b_muon_combRelIso);
   fChain->SetBranchAddress("muon.iSubdet", muon_iSubdet, &b_muon_iSubdet);
   fChain->SetBranchAddress("muon.nMatchedStations", muon_nMatchedStations, &b_muon_nMatchedStations);
   fChain->SetBranchAddress("muon.nLayersWithMmt", muon_nLayersWithMmt, &b_muon_nLayersWithMmt);
   fChain->SetBranchAddress("muon.nValidMuonHits", muon_nValidMuonHits, &b_muon_nValidMuonHits);
   fChain->SetBranchAddress("muon.nValidPixelHits", muon_nValidPixelHits, &b_muon_nValidPixelHits);
   fChain->SetBranchAddress("muon.isGlobalMuon", muon_isGlobalMuon, &b_muon_isGlobalMuon);
   fChain->SetBranchAddress("muon.isPFMuon", muon_isPFMuon, &b_muon_isPFMuon);
   fChain->SetBranchAddress("muon.hasInnerTrack", muon_hasInnerTrack, &b_muon_hasInnerTrack);
   fChain->SetBranchAddress("muon.hasGlobalTrack", muon_hasGlobalTrack, &b_muon_hasGlobalTrack);
   fChain->SetBranchAddress("muon.hasBestTrack", muon_hasBestTrack, &b_muon_hasBestTrack);
   fChain->SetBranchAddress("muon.isLoose", muon_isLoose, &b_muon_isLoose);
   fChain->SetBranchAddress("muon.isTight", muon_isTight, &b_muon_isTight);
   fChain->SetBranchAddress("jet.size", &jet_size, &b_jet_size);
   fChain->SetBranchAddress("jet.pt", jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet.eta", jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet.phi", jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet.px", jet_px, &b_jet_px);
   fChain->SetBranchAddress("jet.py", jet_py, &b_jet_py);
   fChain->SetBranchAddress("jet.pz", jet_pz, &b_jet_pz);
   fChain->SetBranchAddress("jet.energy", jet_energy, &b_jet_energy);
   fChain->SetBranchAddress("jet.jecScale", jet_jecScale, &b_jet_jecScale);
   fChain->SetBranchAddress("jet.chFraction", jet_chFraction, &b_jet_chFraction);
   fChain->SetBranchAddress("jet.nhFraction", jet_nhFraction, &b_jet_nhFraction);
   fChain->SetBranchAddress("jet.ceFraction", jet_ceFraction, &b_jet_ceFraction);
   fChain->SetBranchAddress("jet.neFraction", jet_neFraction, &b_jet_neFraction);
   fChain->SetBranchAddress("jet.iSubdet", jet_iSubdet, &b_jet_iSubdet);
   fChain->SetBranchAddress("jet.nConstituents", jet_nConstituents, &b_jet_nConstituents);
   fChain->SetBranchAddress("jet.nCharged", jet_nCharged, &b_jet_nCharged);
   fChain->SetBranchAddress("jet.isLoose", jet_isLoose, &b_jet_isLoose);
   fChain->SetBranchAddress("vertex.size", &vertex_size, &b_vertex_size);
   fChain->SetBranchAddress("vertex.x", vertex_x, &b_vertex_x);
   fChain->SetBranchAddress("vertex.y", vertex_y, &b_vertex_y);
   fChain->SetBranchAddress("vertex.z", vertex_z, &b_vertex_z);
   fChain->SetBranchAddress("vertex.rho", vertex_rho, &b_vertex_rho);
   fChain->SetBranchAddress("vertex.sumPt2", vertex_sumPt2, &b_vertex_sumPt2);
   fChain->SetBranchAddress("vertex.chi2", vertex_chi2, &b_vertex_chi2);
   fChain->SetBranchAddress("vertex.ndof", vertex_ndof, &b_vertex_ndof);
   fChain->SetBranchAddress("vertex.isGood", vertex_isGood, &b_vertex_isGood);
   fChain->SetBranchAddress("photon.dRGen", photon_dRGen, &b_photon_dRGen);
   fChain->SetBranchAddress("photon.genIso", photon_genIso, &b_photon_genIso);
   fChain->SetBranchAddress("photon.nearestGen", photon_nearestGen, &b_photon_nearestGen);
   fChain->SetBranchAddress("photon.dRJet", photon_dRJet, &b_photon_dRJet);
   fChain->SetBranchAddress("photon.dRNextJet", photon_dRNextJet, &b_photon_dRNextJet);
   fChain->SetBranchAddress("photon.dRPF", photon_dRPF, &b_photon_dRPF);
   fChain->SetBranchAddress("photon.nearestPF", photon_nearestPF, &b_photon_nearestPF);
   fChain->SetBranchAddress("photon.pfIsPU", photon_pfIsPU, &b_photon_pfIsPU);
   fChain->SetBranchAddress("electron.dRGen", electron_dRGen, &b_electron_dRGen);
   fChain->SetBranchAddress("electron.genIso", electron_genIso, &b_electron_genIso);
   fChain->SetBranchAddress("electron.nearestGen", electron_nearestGen, &b_electron_nearestGen);
   fChain->SetBranchAddress("electron.dRJet", electron_dRJet, &b_electron_dRJet);
   fChain->SetBranchAddress("electron.dRNextJet", electron_dRNextJet, &b_electron_dRNextJet);
   fChain->SetBranchAddress("electron.dRPhoton", electron_dRPhoton, &b_electron_dRPhoton);
   fChain->SetBranchAddress("electron.dRNextPhoton", electron_dRNextPhoton, &b_electron_dRNextPhoton);
   fChain->SetBranchAddress("electron.dRPF", electron_dRPF, &b_electron_dRPF);
   fChain->SetBranchAddress("electron.nearestPF", electron_nearestPF, &b_electron_nearestPF);
   fChain->SetBranchAddress("electron.pfIsPU", electron_pfIsPU, &b_electron_pfIsPU);
   fChain->SetBranchAddress("muon.dRGen", muon_dRGen, &b_muon_dRGen);
   fChain->SetBranchAddress("muon.genIso", muon_genIso, &b_muon_genIso);
   fChain->SetBranchAddress("muon.nearestGen", muon_nearestGen, &b_muon_nearestGen);
   fChain->SetBranchAddress("muon.dRJet", muon_dRJet, &b_muon_dRJet);
   fChain->SetBranchAddress("muon.dRNextJet", muon_dRNextJet, &b_muon_dRNextJet);
   fChain->SetBranchAddress("muon.dRPhoton", muon_dRPhoton, &b_muon_dRPhoton);
   fChain->SetBranchAddress("muon.dRNextPhoton", muon_dRNextPhoton, &b_muon_dRNextPhoton);
   fChain->SetBranchAddress("muon.dRPF", muon_dRPF, &b_muon_dRPF);
   fChain->SetBranchAddress("muon.nearestPF", muon_nearestPF, &b_muon_nearestPF);
   fChain->SetBranchAddress("muon.pfIsPU", muon_pfIsPU, &b_muon_pfIsPU);
   fChain->SetBranchAddress("jet.dRGen", jet_dRGen, &b_jet_dRGen);
   fChain->SetBranchAddress("jet.genSumPt", jet_genSumPt, &b_jet_genSumPt);
   fChain->SetBranchAddress("jet.nearestGen", jet_nearestGen, &b_jet_nearestGen);
   Notify();
}

Bool_t allObjects::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void allObjects::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t allObjects::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef allObjects_cxx

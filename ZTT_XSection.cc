////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   Compiling the code:   ./Make.sh ZTT_XSection.cc
//   Running the code:     ./ZTT_XSection.exe OutPut.root   Input.root
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "TreeReader.h"
#include "WeightCalculator.h"
#include <string> //remove space for successful compilation
#include <ostream> //remove space for successful compilation
#include <ctime>

void registerBranches(TTree *tree);
void initBranch(TTree *tree, std::string name, void *add){
  const char *bname = name.c_str();
  tree->SetBranchStatus(bname, 1);
  tree->SetBranchAddress(bname, add);
}

int main(int argc, char** argv) {
  using namespace std;

  std::string out = *(argv + 1);
    
  cout << "OUTPUT: " << out << endl;     //PRINTING THE OUTPUT FILE NAME
  TFile *fout = TFile::Open(out.c_str(), "RECREATE");
    
  std::string input = *(argv + 2);
  cout << " INPUT: " << input << endl;     //PRINTING THE INPUT FILE NAME
  TFile * myFile = TFile::Open(input.c_str());
  TH1F * HistoTot = (TH1F*) myFile->Get("hcount");

  cout<<"nargs: "<<argc<<endl;
  bool splitTauTau = false;
  bool selTauTau = false;
  if(argc>3){
    std::string splitTauTauString = *(argv + 3);
    std::string selTauTauString = *(argv + 4);
    splitTauTau = splitTauTauString == "1";
    selTauTau   = selTauTauString   == "1";
  }

  bool debug = false;
  if(argc>5){
    std::string debugString = *(argv + 5);
    debug = debugString == "1";
  }

  std::cout<<"splitTauTau "<<splitTauTau<<std::endl;
  std::cout<<"selTauTau "<<selTauTau<<std::endl;
  std::cout<<"debug "<<debug<<std::endl;


  //add the histrograms of muon and tau visible mass (both for opposite sign and same sign pair )
  TH1F *    visibleMassOS = new TH1F ("visibleMassOS","visibleMassOS", 30, 0, 300);
  TH1F *    visibleMassSS = new TH1F ("visibleMassSS","visibleMassSS", 30, 0, 300);
  visibleMassOS->SetDefaultSumw2();
  visibleMassSS->SetDefaultSumw2();
    
    
  TTree *Run_Tree = (TTree*) myFile->Get("EventTree");
  cout.setf(ios::fixed, ios::floatfield);

  registerBranches(Run_Tree);

  float MuMass= 0.10565837;
  float eleMass= 0.000511;
  
  //Luminosity
  float LumiWeight = 1;
  if (HistoTot) LumiWeight = weightCalc(HistoTot, input);
  cout << "LumiWeight is " << LumiWeight << "\n";

  //Pileup Reweighting
  TFile * PUData= new TFile("MyDataPileupHistogram2016.root");
  TH1F * HistoPUData= (TH1F *) PUData->Get("pileup");
  HistoPUData->Scale(1.0/HistoPUData->Integral());
  TFile * PUMC= new TFile("mcMoriondPU.root");
  TH1F * HistoPUMC= (TH1F *) PUMC->Get("pileup");
  HistoPUMC->Scale(1.0/HistoPUMC->Integral());

  //cutflow counts
  int nPassTrigger = 0;
  int nPassMuon    = 0;
  int nPassDYVeto  = 0;
  int nPassTau     = 0;
  int nPassMuMET   = 0;
  int nPassBVeto   = 0;

  //keep track of processing time
  std::clock_t start = std::clock();
  double duration;

  // Loop over all events in the TTree.
  auto nentries_wtn = Run_Tree->GetEntries();
  for(auto ievt = 0; ievt < nentries_wtn; ievt++){
    //Print Status
    if((ievt+1)%1000 == 0) {
      duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
      fprintf(stdout, "\r  Processed events: %8d of %8d (%6f events/s)", ievt+1, nentries_wtn, (ievt+1)/duration);
    }
    fflush(stdout);

    Run_Tree->GetEntry(ievt);  // read this event

    //Check if we want to only look at truth tau tau events
    if(splitTauTau){
      int numTau=0;
      for(int igen=0;igen < nMC; igen++){
	if(fabs(mcPID->at(igen)) == 15 && mcMomPID->at(igen)==23) numTau++;
      }
      if( selTauTau && (numTau < 1)) continue; // uncomment this line to get Ztautau contribution of DY
      if(!selTauTau && (numTau > 1)) continue; // uncomment this line to get Zll contribution of DY 
    }

    //
    //Trigger
    //
    bool PassTrigger = (HLTEleMuX >> 19 & 1) == 1; //       else if (name.find("HLT_IsoMu24_v") != string::npos) bitEleMuX = 19;
    if (!PassTrigger) continue;
    nPassTrigger += 1;

    //Compute Pileup Weight
    float PUWeight = 1;
    if (!isData){
      if(debug) std::cout<<"get PUWeight"<<std::endl;
      int puNUmmc=int(puTrue->at(0)*10);
      int puNUmdata=int(puTrue->at(0)*10);
      float PUMC_=HistoPUMC->GetBinContent(puNUmmc+1);
      float PUData_=HistoPUData->GetBinContent(puNUmdata+1);
      PUWeight= PUData_/PUMC_;
    }

    //
    //Loop over muons to find number passing quality cuts. 
    //
    int nSelMuons = 0;
    int mu;
    TLorentzVector Mu4Momentum;
    if(debug) std::cout<<"muon loop 1"<<std::endl;
    for  (int imu=0 ; imu < nMu; imu++){
      //check pt and eta
      if(muPt->at(imu) < 30 || fabs(muEta->at(imu)) > 2.1) continue;

      //compute isolation
      float IsoMu=muPFChIso->at(imu)/muPt->at(imu);
      if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
	IsoMu= ( muPFChIso->at(imu)/muPt->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
      //check isolation
      if(IsoMu > 0.15) continue;

      //Check Muon ID
      if( !(muIDbit->at(imu)>>2 & 1)) continue; // 2 is tight
      nSelMuons += 1;
      mu = imu;
      if(debug) std::cout<<"found muon"<<std::endl;
      Mu4Momentum.SetPtEtaPhiM(muPt->at(mu),muEta->at(mu),muPhi->at(mu),MuMass);
    }
    //Found at least one good muon
    if(nSelMuons == 0) continue;
    nPassMuon += 1;

    //
    // Loose Muon selection for Z->mu mu veto
    //
    bool passDYVeto = true;
    if(debug) std::cout<<"muon loop 2"<<std::endl;
    for(int imu=0 ; imu < nMu; imu++){
      if(imu == mu) continue; //skip the primary muon
      if(muPt->at(imu) < 15 || fabs(muEta->at(imu)) > 2.4) continue;//this one doesnt need to satisfy the muon trigger, so have extended eta range and lower pt threshold

      //compute isolation
      float IsoMu2=muPFChIso->at(imu)/muPt->at(imu);
      if ( (muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu) )  > 0.0)
	IsoMu2= ( muPFChIso->at(imu)/muPt->at(imu) + muPFNeuIso->at(imu) + muPFPhoIso->at(imu) - 0.5* muPFPUIso->at(imu))/muPt->at(imu);
      bool MuIdIso2=((muIDbit->at(imu) >> 0 & 1)  && IsoMu2 < 0.30 && fabs(muD0->at(imu)) < 0.045 && fabs(muDz->at(imu)) < 0.2);
      if(!MuIdIso2) continue;

      //Only worry about the second muon if it is opposite sign of the main muon
      bool OS = muCharge->at(imu) * muCharge->at(mu) < 0;
      if(OS) passDYVeto = false;
    }
    //Only consider events with one muon candidate
    if(passDYVeto != 1) continue;
    nPassDYVeto += 1;

    //
    //Now process Tau's
    //
    int nSelTaus = 0;
    int tau;
    TLorentzVector Tau4Momentum;
    if(debug) std::cout<<"tau loop 1"<<std::endl;
    for  (int itau=0 ; itau < nTau; itau++){
      //check pt and eta
      if(tauPt->at(itau) < 30 || fabs(tauEta->at(itau)) > 2.3) continue;

      //Tau ID
      if(taupfTausDiscriminationByDecayModeFinding->at(itau) < 0.5) continue;
      if(tauByTightMuonRejection3                 ->at(itau) < 0.5) continue;
      if(tauByMVA6LooseElectronRejection          ->at(itau) < 0.5) continue;
      if(tauByTightIsolationMVArun2v1DBoldDMwLT   ->at(itau) < 0.5) continue;
      Tau4Momentum.SetPtEtaPhiM(tauPt->at(itau),tauEta->at(itau),tauPhi->at(itau),tauMass->at(itau));
      if(Tau4Momentum.DeltaR(Mu4Momentum) < 0.5) continue;
      nSelTaus += 1;
      tau = itau;
    } // End of tau loop

    // Only consider events with 1 tau candidate for now
    if(nSelTaus != 1) continue;
    nPassTau += 1;

    //Reject W+Jets
    float MuMetTranverseMass= TMass_F(muPt->at(mu), muPt->at(mu)*cos(muPhi->at(mu)),muPt->at(mu)*sin(muPhi->at(mu)) ,  pfMET, pfMETPhi);      
    if(MuMetTranverseMass > 40) continue;
    nPassMuMET += 1;

    //ttbar veto with bjet veto
    TLorentzVector Jet4Momentum;
    std::vector<float> EveBJetPt;
    EveBJetPt.clear();
    for (int ijet= 0 ; ijet < nJet ; ijet++){
      Jet4Momentum.SetPtEtaPhiE(jetPt->at(ijet),jetEta->at(ijet),jetPhi->at(ijet),jetRawEn->at(ijet));
      if(     jetPt ->at(ijet)  < 20) continue;//jet Pt
      if(fabs(jetEta->at(ijet)) > 2.5) continue;//jet Eta
      if(Jet4Momentum.DeltaR(Tau4Momentum) < 0.5) continue;//Tau Overlap
      if(Jet4Momentum.DeltaR( Mu4Momentum) < 0.5) continue;//Muon Overlap
      if(jetCSV2BJetTags->at(ijet) < 0.8484) continue;//btag
      EveBJetPt.push_back(jetPt->at(ijet));
    }
    //veto events with bjets passing our quality/overlap cuts
    if(EveBJetPt.size()>0) continue;
    nPassBVeto += 1;

    // Check charge of the muon and Taus
    bool  OS = muCharge->at(mu) * tauCharge->at(tau) < 0;
    bool  SS = muCharge->at(mu) * tauCharge->at(tau) > 0;

    // Construct the visible mu+tau system
    TLorentzVector MuTau4Momentum = Mu4Momentum + Tau4Momentum;
    float eventWeight = LumiWeight*PUWeight;

    //Check if there is an OS  muTau pair with dR > 0.5 and TMass(mu.MET) < 40 and then fill the weighted histogram as below:
    if(OS)
      visibleMassOS->Fill(MuTau4Momentum.M(), eventWeight);

    //Check if there is a SS  muTau pair with dR > 0.5 and TMass(mu.MET) < 40 and then fill the weighted histogram as below:
    if(SS)
      visibleMassSS->Fill(MuTau4Momentum.M(), eventWeight);

  } //End Event Loop

  //end of analysis code, close and write histograms/file
  fout->cd();
  visibleMassOS->Write();
  visibleMassSS->Write();
  //visibleMassOSRelaxedTauIso->Write();
  //visibleMassSSRelaxedTauIso->Write();
  fout->Close();

  std::cout << "Done" << std::endl;
  std::cout << "nPassTrigger: " << nPassTrigger << std::endl;
  std::cout << "nPassMuon: " << nPassMuon << std::endl;
  std::cout << "nPassDYVeto: " << nPassDYVeto << std::endl;
  std::cout << "nPassTau: " << nPassTau << std::endl;
  std::cout << "nPassMuMET: " << nPassMuMET << std::endl;
  std::cout << "nPassBVeto: " << nPassBVeto << std::endl;

}



void registerBranches(TTree *tree) {
  tree->SetBranchStatus("*", 0);


  //########################################   General Info
  initBranch(tree, "isData", &isData);
  initBranch(tree,"run", &run);
  initBranch(tree,"lumis", &lumis);
  initBranch(tree,"event", &event);
  initBranch(tree,"genWeight",&genWeight);
  initBranch(tree,"HLTEleMuX", &HLTEleMuX);
  initBranch(tree,"puTrue", &puTrue);
  initBranch(tree,"nVtx",&nVtx);
        
  //########################################   MC Info
  initBranch(tree,"nMC", &nMC);
  initBranch(tree,"mcPID", &mcPID);
  initBranch(tree,"mcStatus", &mcStatus);
  initBranch(tree,"mcPt", &mcPt );
  initBranch(tree,"mcEta", &mcEta );
  initBranch(tree,"mcPhi", &mcPhi );
  initBranch(tree,"mcE", &mcE );
  initBranch(tree,"mcMass", &mcMass );
  initBranch(tree,"mcMomPID", &mcMomPID );
  initBranch(tree,"mcGMomPID", &mcGMomPID );
        
  //########################################   Tau Info
  initBranch(tree,"nTau", &nTau);
  initBranch(tree,"tauPt"  ,&tauPt);
  initBranch(tree,"tauEta"  ,&tauEta);
  initBranch(tree,"tauPhi"  ,&tauPhi);
  initBranch(tree,"tauMass"  ,&tauMass);
  initBranch(tree,"tauCharge"  ,&tauCharge);

  initBranch(tree,"taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding);

  initBranch(tree,"tauByTightMuonRejection3", &tauByTightMuonRejection3);
  initBranch(tree,"tauByLooseMuonRejection3", &tauByLooseMuonRejection3);

  initBranch(tree,"tauByMVA6TightElectronRejection"  ,&tauByMVA6TightElectronRejection);
  initBranch(tree,"tauByMVA6MediumElectronRejection"  ,&tauByMVA6MediumElectronRejection);
  initBranch(tree,"tauByMVA6LooseElectronRejection", &tauByMVA6LooseElectronRejection);

  initBranch(tree,"tauDxy",&tauDxy);
  initBranch(tree,"tauDecayMode",&tauDecayMode);

  initBranch(tree,"tauByLooseIsolationMVArun2v1DBoldDMwLT",&tauByLooseIsolationMVArun2v1DBoldDMwLT);
  initBranch(tree,"tauByVLooseIsolationMVArun2v1DBoldDMwLT",&tauByVLooseIsolationMVArun2v1DBoldDMwLT);
  initBranch(tree,"tauByTightIsolationMVArun2v1DBoldDMwLT",&tauByTightIsolationMVArun2v1DBoldDMwLT);
        
  //########################################   Mu Info
  initBranch(tree,"nMu", &nMu);
  initBranch(tree,"muPt"  ,&muPt);
  initBranch(tree,"muEta"  ,&muEta);
  initBranch(tree,"muPhi"  ,&muPhi);
  initBranch(tree,"muIsoTrk", &muIsoTrk);
  initBranch(tree,"muCharge",&muCharge);
  initBranch(tree,"muIDbit",&muIDbit);//NEW
  initBranch(tree,"muPFChIso", &muPFChIso);
  initBranch(tree,"muPFPhoIso", &muPFPhoIso);
  initBranch(tree,"muPFNeuIso", &muPFNeuIso);
  initBranch(tree,"muPFPUIso", &muPFPUIso);
  initBranch(tree,"muD0",&muD0);
  initBranch(tree,"muDz",&muDz);
        
  //########################################   Ele Info
  initBranch(tree,"nEle", &nEle);
  initBranch(tree,"elePt"  ,&elePt);
  initBranch(tree,"eleEta"  ,&eleEta);
  initBranch(tree,"elePhi"  ,&elePhi);
  initBranch(tree,"elePFChIso", &elePFChIso);
  initBranch(tree,"eleIDMVA", &eleIDMVA);//NEW
  initBranch(tree,"eleCharge",&eleCharge);
  initBranch(tree,"eleSCEta",&eleSCEta);
  initBranch(tree,"elePFChIso", &elePFChIso);
  initBranch(tree,"elePFPhoIso", &elePFPhoIso);
  initBranch(tree,"elePFNeuIso", &elePFNeuIso);
  initBranch(tree,"elePFPUIso", &elePFPUIso);
  initBranch(tree,"eleD0",&eleD0);
  initBranch(tree,"eleDz",&eleDz);
  initBranch(tree,"eleMissHits", &eleMissHits);
  initBranch(tree,"eleConvVeto", &eleConvVeto);
  initBranch(tree,"eleSCEta", &eleSCEta );
        
  //########################################   Jet Info
  initBranch(tree,"nJet",&nJet);
  initBranch(tree,"jetPt",&jetPt);
  initBranch(tree,"jetEta",&jetEta);
  initBranch(tree,"jetPhi",&jetPhi);
  initBranch(tree,"jetCSV2BJetTags",&jetCSV2BJetTags);
  initBranch(tree,"jetPFLooseId",&jetPFLooseId);
  initBranch(tree,"jetPUID",&jetPUID);
  initBranch(tree,"jetRawPt",&jetRawPt);
  initBranch(tree,"jetJECUnc",&jetJECUnc);
  initBranch(tree,"jetRawEn",&jetRawEn);
  initBranch(tree,"jetHadFlvr",&jetHadFlvr);
        
  //########################################   MET Info
  initBranch(tree,"pfMET",&pfMET);
  initBranch(tree,"pfMETPhi",&pfMETPhi);
  initBranch(tree,"metFilters",&metFilters);
  initBranch(tree,"genHT",&genHT);
 

}

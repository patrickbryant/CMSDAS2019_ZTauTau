#ifndef  HISTOGRAMMGR
#define HISTOGRAMMGR

#include <string>
#include <TH1D.h>
#include "TreeReader.h"

std::string XaxisTitle(
    const std::string& title, const std::string& unit );

std::string YaxisTitle(
    const std::string& unit,
    const int nbins, const double xmin, const double xmax  );

///////////////////////////////////////////////////////////////////

class ObjectHisto{
public:
  ObjectHisto( const std::string& cutlevel, const bool domuon );
  ~ObjectHisto();

  void Fill( const double eventweight, const int tau = -1 , const int lep=-1 );
  void Write();
private:
  const std::string cutlevel;
  const bool doMuon;
  std::map<std::string, TH1D> histmap;

  TH1D& Hist( const std::string& );

  void defineHist(
    const std::string& xname,
    const std::string& xtitle,
    const std::string& xunit,
    const unsigned     nbins,
    const double       xmin,
    const double       xmax
  );

  std::string histname( const std::string& );
};

/*** Begin Histogram definition */

ObjectHisto::ObjectHisto( const std::string& cutLevel, const bool domuon ):
 cutlevel( cutLevel ), doMuon(domuon)
{
  defineHist("nJets", "Num. of Jets", "", 10, 0,10);
  defineHist("JetPt", "Jet p_{T}",          "GeV", 50,   0, 500 );
  defineHist("Jet1Pt", "Leading Jet p_{T}", "GeV", 60,   0, 600 );
  defineHist("JetEta", "Jet #eta",          "",    48,-2.4, 2.4 );
  defineHist("Jet1Eta","Leading Jet #eta",  "",    48,-2.4, 2.4 );

  defineHist("MuonPt",  "#mu p_{T}",             "GeV", 30, 0, 300 );
  defineHist("MuonEta", "#mu #eta",              "", 48, -2.4, 2.4 );
  defineHist("ElecPt",  "e p_{T}",              "GeV", 30, 0, 300 );
  defineHist("ElecEta", "e #eta",               "", 48, -2.4, 2.4 );
  defineHist("LepPt",   "Primary lepton p_{T}", "GeV", 48,-2.4,2.4);
  defineHist("LepEta",  "Primary lepton #eta",  "", 48,-2.4,2.4);

  defineHist("MET","Missing Transverse Energy", "GeV", 30,0,300);
  defineHist("TMass", "Lepton transverse mass", "GeV", 50,0,150);
}

void ObjectHisto::Fill( const double eventweight, const int tau, const int iLep )
{
  // Setting up variables
  static const double muMass = 0.10565837;
  static const double elMass = 0.000511  ;
  TLorentzVector lepP4, tauP4;
  if( tau >= 0 )
    tauP4.SetPtEtaPhiM( tauPt->at(tau), tauEta->at(tau), tauPhi->at(tau), tauMass->at(tau) );
  if( doMuon && iLep >= 0 )
    lepP4.SetPtEtaPhiM( muPt->at(iLep), muEta->at(iLep), muPhi->at(iLep), muMass );
  else if ( iLep >= 0 )
    lepP4.SetPtEtaPhiM( elePt->at(iLep), eleEta->at(iLep), elePhi->at(iLep), elMass );
  const double tmass = TMass_F( lepP4.Pt(), lepP4.Px(), lepP4.Py(), pfMET, pfMETPhi );

  // Filling jet histograms
  for( int i = 0 ; i < nJet ; ++i ){
    Hist("JetPt").Fill( jetPt->at(i), eventweight );
    Hist("JetEta").Fill( jetEta->at(i), eventweight );

    if( i == 0 ){
      Hist("Jet1Pt").Fill( jetPt->at(i), eventweight );
    }
  }

  // Muon histograms
  for( int i = 0 ; i < nMu ; ++i ){
    Hist("MuonPt").Fill( muPt->at(i), eventweight );
    Hist("MuonEta").Fill( muEta->at(i), eventweight );
  }

  for( int i = 0 ; i < nEle; ++i ){
    Hist("ElecPt").Fill( elePt->at(i), eventweight );
    Hist("ElecEta").Fill( eleEta->at(i), eventweight );
  }

  if(iLep >=0 ){
    Hist("ElecPt").Fill( lepP4.Pt(), eventweight );
    Hist("ElecEta").Fill( lepP4.Eta(), eventweight );
    Hist("TMass").Fill( tmass, eventweight );
  }

  Hist("MET").Fill(pfMET);
  // TBC
}

/** END Histogram definition */
ObjectHisto::~ObjectHisto() {}

void ObjectHisto::Write()
{
  for( const auto& p : histmap ){
    p.second.Write();
  }
}

TH1D& ObjectHisto::Hist( const std::string& x )
{
  return histmap.at(histname(x));
}

std::string ObjectHisto::histname( const std::string& xname )
{
  return cutlevel + "_" + xname;
}

void ObjectHisto::defineHist(
    const std::string& xname,
    const std::string& xtitle,
    const std::string& xunit,
    const unsigned     nbins,
    const double       xmin,
    const double       xmax )

{
  const std::string name = histname(xname);
  histmap[name] = TH1D( name.c_str(), name.c_str(), nbins, xmin, xmax );
  histmap[name].GetXaxis()->SetTitle( XaxisTitle(xtitle,xunit).c_str() );
  histmap[name].GetYaxis()->SetTitle( YaxisTitle(xunit,nbins,xmax,xmin).c_str() );
}

/////////////////////////////////////////////////////////////

std::string XaxisTitle(
    const std::string& title, const std::string& unit )
{
  std::string ans = title ;
  if( unit != "" ){
    ans += "[" + unit + "]";
  }
  return ans;
}

std::string YaxisTitle(
    const std::string& unit,
    const int nbins, const double xmin, const double xmax  )
{
  char ans[1024];
  const double binwidth = (xmax - xmin)/nbins;

  if( binwidth == 1 ){
    if( unit == "" ){
      return "Events";
    } else {
      sprintf(ans, "Events / %s", unit.c_str() );
    }
  } else {
    sprintf( ans, "Events / %2.2lf %s", binwidth, unit.c_str()  );
  }
  return ans;
}

class EventHisto {
public:
  EventHisto( const std::string& cutlevel, const bool );
  ~EventHisto();

  // Selected Tau, Muon and Electron index (leave -1 is not of this channel)
  void Fill(
      const int tau,
      const int iLep,
      const double evtweight );

  void Write();

private:
  const std::string cutlevel;
  const bool doMuon ;
  std::map<std::string, TH1D> histmap;
  
  TH1D& Hist( const std::string&, const bool, const bool );

  void defineHist(
    const std::string& xname,
    const std::string& xtitle,
    const std::string& xunit,
    const unsigned     nbins,
    const double       xmin,
    const double       xmax
  );

  std::string histname(
    const std::string& name,
    const bool lepiso,
    const bool os
    );
};

/** Begin Editing */

EventHisto::EventHisto( const std::string& x, const bool y ):
cutlevel(x), doMuon(y)
{
  defineHist("visMass", "visible mass", "GeV", 30 , 0 , 300  );
  defineHist("visPt",   "visible p_{T}", "GeV", 30 , 0 , 300  );
  defineHist("visEta",  "visible #eta", "",  24 , -2.4 , 2.4  );

  defineHist("isoMu",    "Muon Isolation", "", 20, 0 , 1 );
  defineHist("LepPt",    "Lepton p_{T}", "GeV", 20, 0 , 200 );
  defineHist("LepEta",   "Lepton #eta", "", 24, -2.4 , 2.4 );
  defineHist("TauPt",    "#tau p_{T}", "GeV", 20, 0 , 200 );
  defineHist("TauEta",   "#tau p_{T}", "GeV", 20, 0 , 200 );
  defineHist("MET",      "MET", "GeV", 20, 0 , 200 );
  defineHist("LepTMass", "Lepton Transverse Mass", "GeV", 20, 0, 150 );
  defineHist("DeltaRLep", "#Delta R(#tau,l)", "", 62,0,6.2);
  defineHist("DeltaRJet", "#Delta R_{min}(#tau,jet)", "",62,0,6.2);
  defineHist("DeltaPhiTMET","#Delta #phi(#tau,MET)", "", 62,-3.2,3.2);
  defineHist("DeltaPhiLepMET", "#Delta #phi(l,MET)", "", 62,-3.2,3.2);
}

void EventHisto::Fill( const int tau, const int iLep, const double evtweight )
{
  static const double muMass = 0.10565837;
  static const double elMass = 0.000511  ;

  TLorentzVector lepP4, tauP4, METP4;
  tauP4.SetPtEtaPhiM( tauPt->at(tau), tauEta->at(tau), tauPhi->at(tau), tauMass->at(tau) );
  bool lepiso = false;
  float IsoMu = 1e6;
  if( doMuon ){
    lepP4.SetPtEtaPhiM( muPt->at(iLep), muEta->at(iLep), muPhi->at(iLep), muMass );
    //compute isolation
    IsoMu  = muPFChIso->at(iLep) / muPt->at(iLep);
    IsoMu += std::max(0., (muPFNeuIso->at(iLep) + muPFPhoIso->at(iLep) - 0.5 * muPFPUIso->at(iLep)) / muPt->at(iLep));
    //check isolation
    lepiso = (IsoMu < 0.15);//tight
  }else{
    lepP4.SetPtEtaPhiM( elePt->at(iLep), eleEta->at(iLep), elePhi->at(iLep), elMass );
    lepiso =
      fabs(eleSCEta->at(iLep)) <= 0.8 && eleIDMVA->at(iLep) > 0.941  ? true :
      fabs(eleSCEta->at(iLep)) >  0.8 && fabs(eleSCEta->at(iLep)) <=  1.5 && eleIDMVA->at(iLep) >   0.899 ? true :
      fabs(eleSCEta->at(iLep)) >=  1.5 && eleIDMVA->at(iLep) >  0.758 ?  true :
      false;
  }
  const TLorentzVector zP4 = tauP4 + lepP4;
  METP4.SetPtEtaPhiE(pfMET,0,pfMETPhi,pfMET);

  //const bool tauiso = tauByTightIsolationMVArun2v1DBoldDMwLT ->at(tau) > 0.5 ;

  const bool os =
    doMuon ?  ( muCharge->at(iLep) * tauCharge->at(tau)  < 0 ) :
              ( eleCharge->at(iLep) * tauCharge->at(tau)  < 0 ) ;

  const double tmass = TMass_F( lepP4.Pt(), lepP4.Px(), lepP4.Py(), pfMET, pfMETPhi );

  // Begin filling
  Hist("isoMu", lepiso,os).Fill(IsoMu, evtweight);
  Hist("visMass",lepiso,os).Fill( zP4.M(), evtweight );
  Hist("visPt"  ,lepiso,os).Fill( zP4.Pt(), evtweight );
  Hist("visEta" ,lepiso,os).Fill( zP4.Eta(), evtweight );

  Hist("LepPt", lepiso, os).Fill( lepP4.Pt(),evtweight );
  Hist("LepEta", lepiso, os).Fill( lepP4.Eta(), evtweight );
  Hist("TauPt",lepiso, os).Fill( tauP4.Pt(), evtweight );
  Hist("TauEta", lepiso, os).Fill( tauP4.Pt(), evtweight );
  Hist("LepTMass", lepiso, os).Fill( tmass, evtweight );

  // angle variables
  Hist("DeltaRLep",lepiso,os).Fill( tauP4.DeltaR(lepP4), evtweight );
  Hist("DeltaPhiTMET",lepiso,os).Fill( tauP4.DeltaPhi(METP4), evtweight );
  Hist("DeltaPhiLepMET",lepiso,os).Fill( lepP4.DeltaPhi(METP4), evtweight );

  double mindeltaR = 1000;
  for( int i = 0 ; i < nJet ; ++i ){
    TLorentzVector jetp4;
    jetp4.SetPtEtaPhiE( jetPt->at(i), jetEta->at(i), jetPhi->at(i), jetEn->at(i));
    mindeltaR = std::min( jetp4.DeltaR( tauP4 ), mindeltaR );
  }
  if( mindeltaR < 1000 ) {
    Hist("DeltaRJet", lepiso, os).Fill(mindeltaR, evtweight );
  }
}

/** Don't Edit After this! */

void EventHisto::defineHist(
    const std::string& xname,
    const std::string& xtitle,
    const std::string& xunit,
    const unsigned     nbins,
    const double       xmin,
    const double       xmax )

{
  for( const auto x : {true,false}){
    for( const auto y : {true,false} ){
      const std::string name = histname(xname,x,y);
      histmap[name] = TH1D( name.c_str(), name.c_str(), nbins, xmin, xmax );
      histmap[name].GetXaxis()->SetTitle( XaxisTitle(xtitle,xunit).c_str() );
      histmap[name].GetYaxis()->SetTitle( YaxisTitle(xunit,nbins,xmin,xmax).c_str() );
    }
  }
}

std::string EventHisto::histname(
  const std::string& x,
  const bool lepiso,
  const bool os
)
{
  const std::string lepisostr = lepiso ? "Iso" : "antiIso" ;
  const std::string osstr     = os ? "OS" : "SS";
  return cutlevel + "_" + x + "_" + lepisostr + "_" + osstr ;
}

TH1D& EventHisto::Hist(
  const std::string& x , const bool lepiso, const bool os )
{
  return histmap.at( histname(x,lepiso,os) );
}


void EventHisto::Write()
{
  for( const auto& hist : histmap ){
    hist.second.Write();
  }
}

EventHisto::~EventHisto() {}


#endif

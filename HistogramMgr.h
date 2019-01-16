#ifndef  HISTOGRAMMGR
#define HISTOGRAMMGR

#include <string> 
#include <TH1D.h>
#include "TreeReader.h"

class EventHisto {
public:
  EventHisto( const std::string& );
  ~EventHisto();

  // Selected Tau, Muon and Electron index (leave -1 is not of this channel)
  void Fill( 
      const int tau, 
      const int muon, 
      const int elec,
      const double evtweight );

  void Write();

private:
  const std::string cutlevel;
  std::map<std::string, TH1D> histmap;
  void defineHist( 
    const std::string& xname,
    const std::string& xtitle,
    const std::string& xunit,
    const unsigned     nbins,
    const double       xmin,
    const double       xmax
  );
};

/** Begin Editing */

EventHisto::EventHisto( const std::string& x ):
  cutlevel(x+"_")
{
  defineHist("visMass", "visible mass", "GeV", 30 , 0 , 300  );
  defineHist("visPt", "visible Pt", "GeV", 30 , 0 , 300  );
  defineHist("visEta", "visible Eta", "",  24 , -2.4 , 2.4  );

  defineHist("LepPt",    "Lepton p_{T}", "GeV", 20, 0 , 200 );
  defineHist("LepEta",   "Lepton #eta", "", 24, -2.4 , 2.4 );
  defineHist("TauPt",    "#tau p_{T}", "GeV", 20, 0 , 200 );
  defineHist("TauEta",   "#tau p_{T}", "GeV", 20, 0 , 200 );
  defineHist("MET",      "MET", "GeV", 20, 0 , 200 );
  defineHist("LepTMass", "Lepton Transverse Mass", "GeV", 25, 0, 50 );
}

void EventHisto::Fill( const int tau, const int mu, const int ele, const double evtweight )
{
  static const double muMass = 0.10565837; 
  static const double elMass = 0.000511  ;
  TLorentzVector tauP4, lepP4; 
  if( tau >= 0 )  tauP4.SetPtEtaPhiM( tauPt->at(tau), tauEta->at(tau), tauPhi->at(tau), tauMass->at(tau) );
  if( mu  >= 0 )  lepP4.SetPtEtaPhiM( muPt->at(mu), muEta->at(mu), muPhi->at(mu), muMass );
  if( ele >= 0 )  lepP4.SetPtEtaPhiM( elePt->at(ele), eleEta->at(ele), elePhi->at(ele), elMass );
  const TLorentzVector zP4 = tauP4 + lepP4; 
  const std::string xs = 
    ( ele >= 0 && (eleCharge->at(ele) * tauCharge->at(tau)) > 0 ) ? "SS" : 
    ( ele >= 0 ) ? "OS"  : 
    ( muCharge->at(mu) * tauCharge->at(tau)  > 0 ) ? "SS" : "OS";

  // Begin filling 
  histmap[cutlevel + "visMass" + xs ].Fill( zP4.M(), evtweight ); 
  histmap[cutlevel + "visPt"  + xs ].Fill( zP4.Pt(), evtweight ); 
  histmap[cutlevel + "visEta" + xs ].Fill( zP4.Eta(), evtweight );

  // TBC
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
  const std::string osname = cutlevel + xname  + "OS";
  const std::string ssname = cutlevel + xname  + "SS";
  histmap[osname] = TH1D( osname.c_str(), osname.c_str(), nbins, xmin, xmax );
  histmap[ssname] = TH1D( ssname.c_str(), ssname.c_str(), nbins, xmin, xmax );
}

void EventHisto::Write()
{
  for( const auto& hist : histmap ){
    hist.second.Write();
  }
}

EventHisto::~EventHisto() {} 

#endif 

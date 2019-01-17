#ifndef  HISTOGRAMMGR
#define HISTOGRAMMGR

#include <string>
#include <TH1D.h>
#include "TreeReader.h"

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
    const bool tauiso,
    const bool os
    );
};

/** Begin Editing */

EventHisto::EventHisto( const std::string& x, const bool y ):
  cutlevel(x), doMuon(y)
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

void EventHisto::Fill( const int tau, const int iLep, const double evtweight )
{
  static const double muMass = 0.10565837;
  static const double elMass = 0.000511  ;

  TLorentzVector tauP4, lepP4;
  tauP4.SetPtEtaPhiM( tauPt->at(tau), tauEta->at(tau), tauPhi->at(tau), tauMass->at(tau) );
  if( doMuon )
    lepP4.SetPtEtaPhiM( muPt->at(iLep), muEta->at(iLep), muPhi->at(iLep), muMass );
  else
    lepP4.SetPtEtaPhiM( elePt->at(iLep), eleEta->at(iLep), elePhi->at(iLep), elMass );
  const TLorentzVector zP4 = tauP4 + lepP4;

  const bool tauiso = tauByTightIsolationMVArun2v1DBoldDMwLT ->at(tau) < 0.5 ;
  const bool os =
    doMuon ?  ( muCharge->at(iLep) * tauCharge->at(tau)  < 0 ) :
              ( eleCharge->at(iLep) * tauCharge->at(tau)  < 0 ) ;


  // Begin filling
  histmap[histname("visMass",tauiso,os)].Fill( zP4.M(), evtweight );
  histmap[histname("visPt"  ,tauiso,os)].Fill( zP4.Pt(), evtweight );
  histmap[histname("visEta" ,tauiso,os)].Fill( zP4.Eta(), evtweight );

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
  const std::string is_os_name = histname(xname,true,true);
  const std::string as_os_name = histname(xname,false,true);
  const std::string is_ss_name = histname(xname,true,false);
  const std::string as_ss_name = histname(xname,false,false);
  histmap[is_os_name] = TH1D( is_os_name.c_str(), is_os_name.c_str(), nbins, xmin, xmax );
  histmap[as_os_name] = TH1D( as_os_name.c_str(), as_os_name.c_str(), nbins, xmin, xmax );
  histmap[is_ss_name] = TH1D( is_ss_name.c_str(), is_ss_name.c_str(), nbins, xmin, xmax );
  histmap[as_ss_name] = TH1D( as_ss_name.c_str(), as_ss_name.c_str(), nbins, xmin, xmax );
}

std::string EventHisto::histname(
  const std::string& x,
  const bool tauiso,
  const bool os
)
{
  const std::string tauisostr = tauiso ? "Iso" : "antiIso" ;
  const std::string osstr     = os ? "OS" : "SS";
  return cutlevel + "_" + x + "_" + tauisostr + "_" + osstr ;
}


void EventHisto::Write()
{
  for( const auto& hist : histmap ){
    hist.second.Write();
  }
}

EventHisto::~EventHisto() {}

#endif

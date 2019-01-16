#ifndef  HISTOGRAMMGR
#defined HISTOGRAMMGR

#include <string> 
#include <TH1D>
#include <TreeReader.h>

enum CutLevel {
  uncut, 
  trigger, 
  nondilepton,
  selectedtau,
  selectedlepton,
  btagveto,
  dummyend
}

class HistoMgr {
public:
  HistoMgr();
  ~HistoMgr();

  // Selected Tau, Muon and Electron index (leave -1 is not of this channel)
  void Fill( const CutLevel, const unsigned int, const unsigned int, const unsigned int);
private:
  std::map<CutLevel, TH1D> histmap;
  void defineHist( 
    const std::string& xname,
    const std::string& xunit,
    const unsigned     nbins,
    const double       xmin,
    const double       xmax
  );
}

/** Begin Editing */

HistoMgr::HistoMgr()
{
  defineHist("Lepton Pt", "GeV", );
}

/** Don't Edit! */



#endif 

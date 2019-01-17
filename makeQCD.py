import ROOT
import copy, sys
import optparse

import math
from array import array
import os
ROOT.gROOT.SetBatch(True)

#sys.path.insert(0, 'XhhResolved/plotting/')
#from plotTools import read_mu_qcd_file

parser = optparse.OptionParser()

parser.add_option('-f', '--flavor',
                  dest='flavor',
                  default="muon"
                  )
parser.add_option('-q', '--qcdFile',
                  dest="qcdFile",
                  default="QCD.root",
                  )

o, a = parser.parse_args()
# if not os.path.isdir(o.outputDir):
#     os.mkdir(o.outputDir)

files={"data": o.flavor+"/data.root",
       "WJets": o.flavor+"/WJetsToLNu.root",
       "TTJets": o.flavor+"/TTJets.root",
       "ZLL": o.flavor+"/DYJetsToLL.root",
       "DY": o.flavor+"/DYJetsToTauTau.root"}


f_data  = ROOT.TFile(files["data"],  "READ")
f_wjets = ROOT.TFile(files["WJets"], "READ")
f_ttbar = ROOT.TFile(files["TTJets"],"READ")
f_zll   = ROOT.TFile(files["ZLL"],   "READ")
f_ztt   = ROOT.TFile(files["DY"],    "READ")

f_qcd   = ROOT.TFile(o.flavor+"/"+o.qcdFile,"RECREATE")



print "Compute SS to OS scale factor using antiIso selection"
#f_qcd   = ROOT.TFile(o.flavor+"/"+o.qcdFile,"READ")
#use a variable that has the full event yield

var="BasicSelection_visMass_antiIso"
SS_h_data   = f_data .Get(var+"_SS")
SS_h_data.SetName("data_"+var+"_SS")
SS_h_wjets  = f_wjets.Get(var+"_SS")
SS_h_wjets.SetName("wjets_"+var+"_SS")
SS_h_ttbar  = f_ttbar.Get(var+"_SS")
SS_h_ttbar.SetName("ttbar_"+var+"_SS")
SS_h_zll  = f_zll.Get(var+"_SS")
SS_h_zll.SetName("zll_"+var+"_SS")
SS_h_ztt  = f_ztt.Get(var+"_SS")
SS_h_ztt.SetName("ztt_"+var+"_SS")

SS_h_qcd   = ROOT.TH1D(SS_h_data)
SS_h_qcd.SetName(var+"_SS")
SS_h_qcd.Add(SS_h_wjets,-1)
SS_h_qcd.Add(SS_h_ttbar,-1)
SS_h_qcd.Add(SS_h_zll,-1)
SS_h_qcd.Add(SS_h_ztt,-1)
SS_h_qcd.Write()

OS_h_data   = f_data .Get(var+"_OS")
OS_h_data.SetName("data_"+var+"_OS")
OS_h_wjets  = f_wjets.Get(var+"_OS")
OS_h_wjets.SetName("wjets_"+var+"_OS")
OS_h_ttbar  = f_ttbar.Get(var+"_OS")
OS_h_ttbar.SetName("ttbar_"+var+"_OS")
OS_h_zll  = f_zll.Get(var+"_OS")
OS_h_zll.SetName("zll_"+var+"_OS")
OS_h_ztt  = f_ztt.Get(var+"_OS")
OS_h_ztt.SetName("ztt_"+var+"_OS")

OS_h_qcd   = ROOT.TH1D(OS_h_data)
OS_h_qcd.SetName(var+"_OS")
OS_h_qcd.Add(OS_h_wjets,-1)
OS_h_qcd.Add(OS_h_ttbar,-1)
OS_h_qcd.Add(OS_h_zll,-1)
OS_h_qcd.Add(OS_h_ztt,-1)
OS_h_qcd.Write()

e_SS = ROOT.Double(0)
e_OS = ROOT.Double(0)
n_SS = SS_h_qcd.IntegralAndError(0, SS_h_qcd.GetSize()-1, e_SS)
n_OS = OS_h_qcd.IntegralAndError(0, OS_h_qcd.GetSize()-1, e_OS)

print "SS =",n_SS,"+/-",e_SS
print "OS =",n_OS,"+/-",e_OS

r = n_OS/n_SS
e = ((e_OS*1.0/n_SS)**2 + (e_SS*n_OS/n_SS**2)**2)**0.5
print "QCD_SS_to_OS_SF =",r,"+/-",e



def subtract():
    # for dName in inFile.GetListOfKeys():
    #     #if "TwoTag" not in dName.GetName(): continue
    #     print dName,dName.GetClassName()
    #     thisDirName = dName.GetName()
    #     dataDir  = inFile.Get(thisDirName)
    #     f_qcd.mkdir(thisDirName)
    #     f_qcd.cd(thisDirName)
    for histKey in f_data.GetListOfKeys():
        histName = histKey.GetName()
        if var in histName: continue #already stored this
        print histName
        h_data   = f_data .Get(histName)
        h_data.SetName("data_"+histName)
        h_wjets  = f_wjets.Get(histName)
        h_wjets.SetName("wjets_"+histName)
        h_ttbar  = f_ttbar.Get(histName)
        h_ttbar.SetName("ttbar_"+histName)
        h_zll  = f_zll.Get(histName)
        h_zll.SetName("zll_"+histName)
        h_ztt  = f_ztt.Get(histName)
        h_ztt.SetName("ztt_"+histName)

        h_qcd   = ROOT.TH1D(h_data)
        h_qcd.SetName(histName)
        h_qcd.Add(h_wjets,-1)
        h_qcd.Add(h_ttbar,-1)
        h_qcd.Add(h_zll,-1)
        h_qcd.Add(h_ztt,-1)
        h_qcd.Scale(r)
        h_qcd.Write()

print "Make QCD hists by subtracting SS MC from SS Data and then scaling by QCD_SS_to_OS_SF =",r,"+/-",e
print " data:",f_data
print "  qcd:",f_qcd
subtract()
f_qcd.Close()


# print "Compute SS to OS scale factor using antiIso selection"
# f_qcd   = ROOT.TFile(o.flavor+"/"+o.qcdFile,"READ")
# #use a variable that has the full event yield

# var="BasicSelection_visMass_antiIso"
# h_SS = f_qcd.Get(var+"_SS")
# h_OS = f_qcd.Get(var+"_OS")

# e_SS = ROOT.Double(0)
# e_OS = ROOT.Double(0)
# n_SS = h_SS.IntegralAndError(0, h_SS.GetSize()-1, e_SS)
# n_OS = h_OS.IntegralAndError(0, h_OS.GetSize()-1, e_OS)

# print "SS =",n_SS,"+/-",e_SS
# print "OS =",n_OS,"+/-",e_OS

# r = n_OS/n_SS
# e = ((e_OS*1.0/n_SS)**2 + (e_SS*n_OS/n_SS**2)**2)**0.5
# print "QCD_SS_to_OS_SF =",r,"+/-",e

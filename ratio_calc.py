# Code written by:  zaixing.mao@cern.ch && edward.laird@cern.ch from Brown U.
#!/usr/bin/env python
import ROOT as r
import numpy
from sys import argv, exit, stdout, stderr
import math

dirName = argv[1]


r.gStyle.SetOptStat(0)
r.gROOT.SetBatch(True)  # to suppress canvas pop-outs


f = r.TFile("prefit.root","recreate")

def getBins(hist, mass_low, mass_high):
    bin_low = -1
    bin_high = -1
    for i in range(hist.GetNbinsX()):
        if hist.GetBinCenter(i+1) >= mass_low and bin_low == -1:
            bin_low = i+1
        if hist.GetBinCenter(i+1) >= mass_high and bin_high == -1:
            bin_high = i
        if bin_low != -1 and bin_high != -1:
            return bin_low, bin_high
    return bin_low, bin_high


def ratio_calculator(fileList = [], mass_low = 25, mass_high = 125, nbins = 50, variableName = "visibleMass", binsSetting = [30, 0, 150]):

   
    IsoOrAnti = ["Iso_","antiIso_"]
    WJets_MC_IA = 0.0 
    WJets_Data_IA = 0.0

    #loop over all the samples
    for iFileName, iFileLocation in fileList:
        isData = False
        if iFileName == 'data':
            isData = True
        isWJet = False
        if iFileName == 'WJets':
            isWJet = True

        print iFileName

        ifile = r.TFile(iFileLocation)

        weight = -1.0
        tauWeight = 0.9
        
        if isData or isWJet:
            weight = 1.0
            tauWeight = 1.0

        for i,iA in enumerate(IsoOrAnti):
           osName = "%sOS" %(variableName+iA)
           ssName = "%sSS" %(variableName+iA)

           lowBin, highBin = getBins(ifile.Get(osName), mass_low, mass_high)

           if isWJet:
               WJets_MC_IA += weight*ifile.Get(osName).Integral(lowBin, highBin)
           else:
               WJets_Data_IA += weight*ifile.Get(osName).Integral(lowBin, highBin)
     
    
    ratio = WJets_Data_IA / WJets_MC_IA
    print("WJet_Data to WJet_MC ratio: "+str(ratio))
    

   

   

fileList = [('DY', '%s/DYJetsToTauTau.root' %dirName),
            ('ZLL', '%s/DYJetsToLL.root ' %dirName),
            ('TTJets', '%s/TTJets.root' %dirName),
            ('WJets', '%s/WJetsToLNu.root' %dirName),
#            ('Diboson', '%s/WZ.root' %dirName),
#            ('Diboson', '%s/ZZ.root' %dirName),
            ('data', '%s/data.root' %dirName),
            ]

variableName = "BeforeTMass_LepTMass_"
bining = []
ratio_calculator(fileList, 80, 160, 30, variableName, [30, 0, 300])

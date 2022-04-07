import ROOT
import os
from MTopCorrelations.Tools.helpers import getObjFromFile

maps_el = {
    2016: {
        'tight' : "egm_preUL/2016/tight/egammaEffi.txt_EGM2D.root",
        'VL'    : "egm_preUL/2016/VL/egammaEffi.txt_EGM2D.root",
    },
    2017: {
        'tight' : "egm_preUL/2017/tight/egammaEffi.txt_EGM2D.root",
        'VL'    : "egm_preUL/2017/VL/egammaEffi.txt_EGM2D.root",
    },
    2018: {
        'tight' : "egm_preUL/2018/tight/egammaEffi.txt_EGM2D.root",
        'VL'    : "egm_preUL/2018/VL/egammaEffi.txt_EGM2D.root",
    }
}

maps_mu = {
    2016: {
        'tight' : "muonNoPs/2016/muonTOPLeptonMVATight2016.root",
        'VL'    : "muonNoPs/2016/muonTOPLeptonMVAVLoose2016.root",
    },
    2017: {
        'tight' : "muonNoPs/2017/muonTOPLeptonMVATight2017.root",
        'VL'    : "muonNoPs/2017/muonTOPLeptonMVAVLoose2017.root",
    },
    2018: {
        'tight' : "muonNoPs/2018/muonTOPLeptonMVATight2018.root",
        'VL'    : "muonNoPs/2018/muonTOPLeptonMVAVLoose2018.root",
    }
}    
    
class leptonSF_topMVA:
    def __init__(self, year, ID = None):
        self.dataDir = "$CMSSW_BASE/src/tWZ/Tools/data/mvaTOP/"
        self.year = year
        filename_el = self.dataDir+maps_el[year][ID]
        filename_mu = self.dataDir+maps_mu[year][ID]
        
        self.SFmaps = {
            "elec" : {
                "SF" :  getObjFromFile(self.dataDir+maps_el[year][ID],"EGamma_SF2D"),
                "sys":  getObjFromFile(self.dataDir+maps_el[year][ID],"sys"),
                "stat": getObjFromFile(self.dataDir+maps_el[year][ID],"stat"),
            },
            "muon" : {
                "SF" :  getObjFromFile(self.dataDir+maps_mu[year][ID],"SF"),
                "sys":  getObjFromFile(self.dataDir+maps_mu[year][ID],"SFTotUnc"),
                "stat": getObjFromFile(self.dataDir+maps_mu[year][ID],"SFTotStat"),                
            }
        }

        
    def getSF(self, pdgId, pt, eta, unc='sys', sigma=0):
        uncert = "sys"
        if unc == "stat":
            uncert = "stat"
        lepton = None
        if abs(pdgId)==11:
            lepton = "elec"
            if eta > 2.5:
                eta = 2.49 
            if eta < -2.5:
                eta = -2.49
            if pt > 200:
                pt = 199

        elif abs(pdgId)==13:
            lepton = "muon"
            if eta > 2.4:
                eta = 2.39 
            if pt > 120:
                pt = 119 
        else: 
          raise Exception("Lepton SF for PdgId %i not known"%pdgId)
          
        etabin = self.SFmaps[lepton]["SF"].GetXaxis().FindBin(eta)
        ptbin  = self.SFmaps[lepton]["SF"].GetYaxis().FindBin(pt)
        
        SF = self.SFmaps[lepton]["SF"].GetBinContent(etabin, ptbin)
        err = self.SFmaps[lepton][uncert].GetBinContent(etabin, ptbin)
        
        
        return SF+sigma*err

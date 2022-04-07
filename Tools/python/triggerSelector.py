import logging
logger = logging.getLogger(__name__)

class triggerSelector:
    # ttH multilepton AN Run II
    # http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2019_111_v6.pdf

    @staticmethod
    def getTriggerList( sample ):
        return [t.GetName() for t in sample.chain.GetListOfBranches() if t.GetName().startswith("HLT_")]

    def __init__(self, year):
        if year == 2016:
            self.mmm    = ["HLT_TripleMu_12_10_5"]
            self.mm     = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"]
            self.m      = ['HLT_IsoMu24','HLT_IsoTkMu24']
            self.eee    = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"]
            self.ee     = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"]
            self.e      = ["HLT_Ele27_WPTight_Gsf"]
            self.em     = ["HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"]
            self.eem    = ["HLT_Mu8_DiEle12_CaloIdL_TrackIdL"]
            self.emm    = ["HLT_DiMu9_Ele9_CaloIdL_TrackIdL"]

        elif year == 2017:
            self.mmm    = ["HLT_TripleMu_10_5_5_DZ", "HLT_TripleMu_12_10_5"]
            self.mm     = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8"]
            #self.m      = ["HLT_IsoMu27", "HLT_IsoMu30"] #TOP-18-009
            self.m      = [ "HLT_IsoMu24", "HLT_IsoMu27", "HLT_IsoMu30"] #Mu24 from ttH multilepton AN 2019/111 v6 
            self.eee    = ["HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"]
            self.ee     = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"]
            self.e      = ["HLT_Ele32_WPTight_Gsf", "HLT_Ele35_WPTight_Gsf", "HLT_Ele38_WPTight_Gsf", "HLT_Ele40_WPTight_Gsf"]
            self.em     = ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"]
            self.eem    = ["HLT_Mu8_DiEle12_CaloIdL_TrackIdL", "HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ"]
            self.emm    = ["HLT_DiMu9_Ele9_CaloIdL_TrackIdL", "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ"]

        elif year == 2018: #from ttH multilepton AN 2019/111 v6. The commented triggers are from 2017 and are available in 2018 as well.
            self.mmm    = [ "HLT_TripleMu_12_10_5" ]
            self.mm     = [ "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"] # HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8
            self.m      = [ "HLT_IsoMu24", "HLT_IsoMu27", "HLT_IsoMu30"] 
            self.eee    = [ "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"]
            self.ee     = [ "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL" ] #HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
            self.e      = [ "HLT_Ele32_WPTight_Gsf" ] # higher threshold not needed of Ele32 was unprescaled in 17+18
            self.em     = [ "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ" ] #"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"
            self.eem    = [ "HLT_Mu8_DiEle12_CaloIdL_TrackIdL"] # HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ
            self.emm    = [ "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ" ] # HLT_DiMu9_Ele9_CaloIdL_TrackId

        else:
            raise NotImplementedError("Trigger selection %r not implemented"%year)

        # define an arbitrary hierarchy
        if year == 2016 or year == 2017:
            self.PDHierarchy = [ "DoubleMuon", "DoubleEG", "MuonEG", "SingleMuon", "SingleElectron" ]
        else:
            # DoubleEG and SingleElectron PDs are merged into EGamma. No change necessary for MC though.
            self.PDHierarchy = [ "DoubleMuon", "EGamma", "MuonEG", "SingleMuon" ]

    def __getVeto(self, cutString):
        return "!%s"%cutString

    def getSelection(self, PD, triggerList = None):

        # reduce the set of triggers in case 
        for lst in [ "mmm", "mm", "m", "eee", "ee", "e", "em", "eem", "emm" ]:
            res = []
            for trig in getattr(self, lst):
                if type(triggerList)==type([]) and not trig in triggerList:
                    logger.warning( "Removing trigger %s", trig )
                else:
                    res.append(trig)
            setattr( self, 'red_'+lst, res )

        # define which triggers should be used for which dataset
        self.DoubleMuon     = "(%s)"%"||".join(self.red_mmm + self.red_mm)
        self.DoubleEG       = "(%s)"%"||".join(self.red_eee + self.red_ee)
        self.MuonEG         = "(%s)"%"||".join(self.red_em + self.red_eem + self.red_emm)
        self.SingleMuon     = "(%s)"%"||".join(self.red_m)
        self.SingleElectron = "(%s)"%"||".join(self.red_e)
        self.EGamma         = "(%s)"%"||".join(self.red_ee+self.red_e)

        found = False
        cutString = ""
        if PD == "MC":
            return "(%s)"%"||".join([self.DoubleMuon, self.DoubleEG, self.MuonEG, self.SingleMuon, self.SingleElectron])
        else:
            for x in reversed(self.PDHierarchy):
                if found:
                    cutString += "&&%s"%self.__getVeto(getattr(self,x))
                if x in PD:# == getattr(self, PD):
                    found = True
                    cutString = getattr(self, x)

            return "(%s)"%cutString

import copy, os, sys
from RootTools.core.Sample import Sample
import ROOT

# Logging
import logging
logger = logging.getLogger(__name__)

from MTopCorrelations.samples.color import color

# Data directory
try:
    directory_ = sys.modules['__main__'].directory_
except:
    import MTopCorrelations.samples.UL_nanoAODv9_locations as locations
    directory_ = locations.mc_UL2018

logger.info("Loading MC samples from directory %s", directory_)

def make_dirs( dirs ):
    return [ os.path.join( directory_, dir_ ) for dir_ in dirs ]

dirs = {}

sample_path_3 = '/groups/hephy/cms/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'
TTbar_3 = Sample.fromDirectory(name="TTbar_172p5", treeName="Events", isData=False, color=color.TTbar, texName="t#bar{t}", directory=sample_path_3)

sample_path_1 = '/groups/hephy/cms/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_mtop169p5_TuneCP5_13TeV-powheg-pythia8'
TTbar_1 = Sample.fromDirectory(name="TTbar_169p5", treeName="Events", isData=False, color=ROOT.kGreen, texName="t#bar{t}", directory=sample_path_1)

sample_path_2 = '/groups/hephy/cms/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8'
TTbar_2 = Sample.fromDirectory(name="TTbar_171p5", treeName="Events", isData=False, color=ROOT.kBlue, texName="t#bar{t}", directory=sample_path_2)

sample_path_4 = '/groups/hephy/cms/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8'
TTbar_4 = Sample.fromDirectory(name="TTbar_173p5", treeName="Events", isData=False, color=ROOT.kYellow, texName="t#bar{t}", directory=sample_path_4)

sample_path_5 = '/groups/hephy/cms/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_mtop175p5_TuneCP5_13TeV-powheg-pythia8'
TTbar_5 = Sample.fromDirectory(name="TTbar_175p5", treeName="Events", isData=False, color=ROOT.kRed, texName="t#bar{t}", directory=sample_path_5)

sample_path_c = '/groups/hephy/cms/dennis.schwarz/MTopCorrelations/QCD_Pt-1000_MuEnrichedPt5_TuneCP5_13TeV-pythia8'
QCD_0 = Sample.fromDirectory(name="QCD_sample", treeName="Events", isData=False, color=ROOT.kYellow, texName="t#bar{t}", directory=sample_path_c)

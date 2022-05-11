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

sample_path_2 = '/scratch-cbe/users/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'
TTbar_2 = Sample.fromDirectory(name="TTbar_172p5", treeName="Events", isData=False, color=color.TTbar, texName="t#bar{t}", directory=sample_path_2)

sample_path_1 = '/scratch-cbe/users/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_mtop171p5_TuneCP5_13TeV-powheg-pythia8'
TTbar_1 = Sample.fromDirectory(name="TTbar_171p5", treeName="Events", isData=False, color=ROOT.kBlue, texName="t#bar{t}", directory=sample_path_1)

sample_path_3 = '/scratch-cbe/users/dennis.schwarz/MTopCorrelations/TTToSemiLeptonic_mtop173p5_TuneCP5_13TeV-powheg-pythia8'
TTbar_3 = Sample.fromDirectory(name="TTbar_173p5", treeName="Events", isData=False, color=ROOT.kYellow, texName="t#bar{t}", directory=sample_path_3)

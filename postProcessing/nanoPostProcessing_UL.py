#!/usr/bin/env python

# standard imports
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import sys
import os
import subprocess
import shutil
import uuid

from array import array
from operator import mul
from math import sqrt, atan2, sin, cos

# RootTools
from RootTools.core.standard import *

# User specific
import tWZ.Tools.user as user

# Tools for systematics
from tWZ.Tools.helpers             import closestOSDLMassToMZ, deltaR, deltaPhi, bestDRMatchInCollection, nonEmptyFile, getSortedZCandidates, cosThetaStar, m3, getMinDLMass
from tWZ.Tools.objectSelection_UL     import getMuons, getElectrons, muonSelector, eleSelector, getGoodMuons, getGoodElectrons, isBJet, getGenPartsAll, getJets, genLepFromZ, mvaTopWP
from tWZ.Tools.objectSelection_UL     import getGenZs

from tWZ.Tools.overlapRemovalTTG   import hasMesonMother, getParentIds
from tWZ.Tools.triggerEfficiency   import triggerEfficiency
from tWZ.Tools.leptonSF            import leptonSF as leptonSF_
from tWZ.Tools.mcTools import pdgToName, GenSearch, B_mesons, D_mesons, B_mesons_abs, D_mesons_abs
genSearch = GenSearch()

from tWZ.Tools.metFiltersUL               import getFilterCut
from Analysis.Tools.puProfileDirDB           import puProfile
from Analysis.Tools.LeptonTrackingEfficiency import LeptonTrackingEfficiency
from Analysis.Tools.helpers                  import checkRootFile, deepCheckRootFile, deepCheckWeight, dRCleaning
from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons
# central configuration
targetLumi = 1000 #pb-1 Which lumi to normalize to

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for cmgPostProcessing")

    argParser.add_argument('--logLevel',    action='store',         nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
    argParser.add_argument('--samples',     action='store',         nargs='*',  type=str, default=['TTZToLLNuNu_ext'],                  help="List of samples to be post-processed, given as CMG component name" )
    argParser.add_argument('--eventsPerJob',action='store',         nargs='?',  type=int, default=30000000,                             help="Maximum number of events per job (Approximate!)." )
    argParser.add_argument('--interpolationOrder',   action='store',nargs='?',  type=int, default=2,                                    help="Maximum number of simultaneous jobs." )
    argParser.add_argument('--nJobs',       action='store',         nargs='?',  type=int, default=2,                                    help="EFT interpolation order" )
    argParser.add_argument('--job',         action='store',                     type=int, default=0,                                    help="Run only jobs i" )
    argParser.add_argument('--minNJobs',    action='store',         nargs='?',  type=int, default=1,                                    help="Minimum number of simultaneous jobs." )
    argParser.add_argument('--targetDir',   action='store',         nargs='?',  type=str, default=user.postprocessing_output_directory, help="Name of the directory the post-processed files will be saved" )
    argParser.add_argument('--processingEra', action='store',       nargs='?',  type=str, default='postProcessed_80X_v22',              help="Name of the processing era" )
    argParser.add_argument('--skim',        action='store',         nargs='?',  type=str, default='trilep',                             help="Skim conditions to be applied for post-processing" )
    argParser.add_argument('--LHEHTCut',    action='store',         nargs='?',  type=int, default=-1,                                   help="LHE cut." )
    argParser.add_argument('--year',        action='store',                     type=str,                                               help="Which year?" )
    argParser.add_argument('--overwriteJEC',action='store',                               default=None,                                 help="Overwrite JEC?" )
    argParser.add_argument('--overwrite',   action='store_true',                                                                        help="Overwrite existing output files, bool flag set to True  if used" )
    argParser.add_argument('--small',       action='store_true',                                                                        help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
    argParser.add_argument('--flagTTBar',   action='store_true',                                                                        help="Is ttbar?" )
    argParser.add_argument('--doCRReweighting',             action='store_true',                                                        help="color reconnection reweighting?")
    argParser.add_argument('--triggerSelection',            action='store_true',                                                        help="Trigger selection?" )
    #argParser.add_argument('--skipGenLepMatching',          action='store_true',                                                        help="skip matched genleps??" )
    argParser.add_argument('--checkTTGJetsOverlap',         action='store_true',                                                        help="Keep TTGJetsEventType which can be used to clean TTG events from TTJets samples" )
    argParser.add_argument('--skipSystematicVariations',    action='store_true',                                                        help="Don't calulcate BTag, JES and JER variations.")
    argParser.add_argument('--noTopPtReweighting',          action='store_true',                                                        help="Skip top pt reweighting.")
    argParser.add_argument('--forceProxy',                  action='store_true',                                                        help="Don't check certificate")
    argParser.add_argument('--skipNanoTools',               action='store_true',                                                        help="Skipt the nanoAOD tools step for computing JEC/JER/MET etc uncertainties")
    argParser.add_argument('--keepNanoAOD',                 action='store_true',                                                        help="Keep nanoAOD output?")
    argParser.add_argument('--reuseNanoAOD',                action='store_true',                                                        help="Keep nanoAOD output?")
    argParser.add_argument('--reapplyJECS',                 action='store_true',                                                        help="Reapply JECs to data?")
    argParser.add_argument('--reduceSizeBy',                action='store',     type=int,                                               help="Reduce the size of the sample by a factor of...")
    argParser.add_argument('--event',                       action='store',     type=int, default=-1,                                   help="Just process event no")

    return argParser

options = get_parser().parse_args()

if options.year not in [ 'UL2016', 'UL2016_preVFP', 'UL2017', 'UL2018' ]:
    raise Exception("Year %s not known"%year)
yearint = 2018 
if "2016" in options.year: 
    yearint = 2016 
elif "2017" in options.year:
    yearint = 2017 
elif "2018" in options.year:
    yearint = 2018

# Logging
import tWZ.Tools.logger as _logger
logFile = '/tmp/%s_%s_%s_njob%s.txt'%(options.skim, '_'.join(options.samples), os.environ['USER'], str(0 if options.nJobs==1 else options.job))
logger  = _logger.get_logger(options.logLevel, logFile = logFile)

#import Analysis.Tools.logger as _logger_an
#logger_an = _logger_an.get_logger(options.logLevel, logFile = logFile )

import RootTools.core.logger as _logger_rt
logger_rt = _logger_rt.get_logger(options.logLevel, logFile = logFile )

def fill_vector_collection( event, collection_name, collection_varnames, objects, maxN = 100):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects[:maxN]):
        for var in collection_varnames:
            if var in obj.keys():
                if type(obj[var]) == type("string"):
                    obj[var] = int(ord(obj[var]))
                if type(obj[var]) == type(True):
                    obj[var] = int(obj[var])
                getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

# Flags
isDiLep         = options.skim.lower().count('dilep') > 0
isTriLep        = options.skim.lower().count('trilep') > 0
isSingleLep     = options.skim.lower().count('singlelep') > 0
isSmall         = options.skim.lower().count('small') > 0
isInclusive     = options.skim.lower().count('inclusive') > 0

# Skim condition
skimConds = []

if options.event > 0:
    skimConds.append( "event == %s"%options.event )

if isDiLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.4)>=2" )
elif isTriLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)&&Electron_pfRelIso03_all<0.4) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5&&Muon_pfRelIso03_all<0.4)>=2 && Sum$(Electron_pt>10&&abs(Electron_eta)<2.5)+Sum$(Muon_pt>10&&abs(Muon_eta)<2.5)>=3" )
elif isSingleLep:
    skimConds.append( "Sum$(Electron_pt>20&&abs(Electron_eta)<2.5) + Sum$(Muon_pt>20&&abs(Muon_eta)<2.5)>=1" )

if isInclusive:
    skimConds.append('(1)')
    isSingleLep = True #otherwise no lepton variables?!
    isDiLep     = True
    isTriLep    = True

################################################################################
#Samples: Load samples
maxN = 1 if options.small else None
if options.small:
    options.job = 0
    options.nJobs = 10000 # set high to just run over 1 input file

if options.year == "UL2016":
    from Samples.nanoAOD.UL16_nanoAODv9         import allSamples as mcSamples
    # allSamples = mcSamples + dataSamples
    allSamples = mcSamples
elif options.year == "UL2016_preVFP": #AODAPV samples
    from Samples.nanoAOD.UL16_nanoAODAPVv9         import allSamples as mcSamples
    # allSamples = mcSamples + dataSamples
    allSamples = mcSamples
elif options.year == "UL2017":
    from Samples.nanoAOD.UL17_nanoAODv9           import allSamples as mcSamples
    # allSamples = mcSamples + dataSamples
    allSamples = mcSamples
elif options.year == "UL2018":
    from Samples.nanoAOD.UL18_nanoAODv9           import allSamples as mcSamples
    # allSamples = mcSamples + dataSamples
    allSamples = mcSamples

samples = []
for selected in options.samples:
    for sample in allSamples:
        if selected == sample.name:
            samples.append(sample)

if len(samples)==0:
    logger.info( "No samples found. Was looking for %s. Exiting" % options.samples )
    sys.exit(-1)

isData = False not in [s.isData for s in samples]
isMC   =  True not in [s.isData for s in samples]

# Check that all samples which are concatenated have the same x-section.
assert isData or len(set([s.xSection for s in samples]))==1, "Not all samples have the same xSection: %s !"%(",".join([s.name for s in samples]))
assert isMC or len(samples)==1, "Don't concatenate data samples"

xSection = samples[0].xSection if isMC else None

# apply MET filter
skimConds.append( getFilterCut(options.year, isData=isData, ignoreJSON=True, skipWeight=True, skipECalFilter=True) )

triggerEff          = triggerEfficiency(yearint)

################################################################################
#Samples: combine if more than one
if len(samples)>1:
    sample_name =  samples[0].name+"_comb"
    logger.info( "Combining samples %s to %s.", ",".join(s.name for s in samples), sample_name )
    sample      = Sample.combine(sample_name, samples, maxN = maxN)
    sampleForPU = Sample.combine(sample_name, samples, maxN = -1)
elif len(samples)==1:
    sample      = samples[0]
    sampleForPU = samples[0]
else:
    raise ValueError( "Need at least one sample. Got %r",samples )

################################################################################
# Reduce samples?
if options.reduceSizeBy > 1:
    logger.info("Sample size will be reduced by a factor of %s", options.reduceSizeBy)
    logger.info("Recalculating the normalization of the sample. Before: %s", sample.normalization)
    if isData:
        NotImplementedError ( "Data samples shouldn't be reduced in size!!" )
    sample.reduceFiles( factor = options.reduceSizeBy )
    # recompute the normalization
    sample.clear()
    sample.name += "_redBy%s"%options.reduceSizeBy
    sample.normalization = sample.getYieldFromDraw(weightString="genWeight")['val']
    logger.info("New normalization: %s", sample.normalization)

################################################################################
options.skim = options.skim + '_small' if options.small else options.skim

################################################################################
# LHE cut (DY samples)
if options.LHEHTCut>0:
    sample.name+="_lheHT"+str(options.LHEHTCut)
    logger.info( "Adding upper LHE cut at %f", options.LHEHTCut )
    skimConds.append( "LHE_HTIncoming<%f"%options.LHEHTCut )

################################################################################
# final output directory
storage_directory = os.path.join( options.targetDir, options.processingEra, options.year, options.skim, sample.name )
try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(storage_directory)
    logger.info( "Created output directory %s.", storage_directory )
except:
    pass

################################################################################
# EFT reweighting
addReweights = False
if hasattr( sample, "reweight_pkl" ):
    addReweights = True
    from Analysis.Tools.WeightInfo                   import WeightInfo
    from Analysis.Tools.HyperPoly                    import HyperPoly
    weightInfo = WeightInfo( sample.reweight_pkl )
    ref_point_coordinates = [ weightInfo.ref_point_coordinates[var] for var in weightInfo.variables ]

    def interpret_weight(weight_id):
        str_s = weight_id.split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res

    hyperPoly       = HyperPoly( options.interpolationOrder )

    weightInfo_data = list(weightInfo.data.iteritems())
    weightInfo_data.sort( key = lambda w: w[1] )
    basepoint_coordinates = map( lambda d: [d[v] for v in weightInfo.variables] , map( lambda w: interpret_weight(w[0]), weightInfo_data) )

    # find the ref point position
    try:
        ref_point_index = basepoint_coordinates.index(ref_point_coordinates)
    except ValueError:
        ref_point_index = None

    hyperPoly.initialize( basepoint_coordinates, ref_point_coordinates )

    logger.info("Adding reweights. Expect to read %i base point weights.", weightInfo.nid)

################################################################################
## sort the list of files?
len_orig = len(sample.files)
sample = sample.split( n=options.nJobs, nSub=options.job)
logger.info( "fileBasedSplitting: Run over %i/%i files for job %i/%i."%(len(sample.files), len_orig, options.job, options.nJobs))
logger.debug("fileBasedSplitting: Files to be run over:\n%s", "\n".join(sample.files) )

################################################################################
# systematic variations
addSystematicVariations = (not isData) and (not options.skipSystematicVariations)

################################################################################
# TODO FIX BTAG REWEIGHTING!!!!

# B tagging SF
# CHECK FOR UL FILES
b_tagger = "DeepJet"
from Analysis.Tools.BTagEfficiency import BTagEfficiency
btagEff = BTagEfficiency( fastSim = False, year=yearint, tagger=b_tagger )

################################################################################
# tmp_output_directory
tmp_output_directory  = os.path.join( user.postprocessing_tmp_directory, "%s_%s_%s_%s_%s"%(options.processingEra, options.year, options.skim, sample.name, str(uuid.uuid3(uuid.NAMESPACE_OID, sample.name))))
if os.path.exists(tmp_output_directory) and options.overwrite:
    if options.nJobs > 1:
        logger.warning( "NOT removing directory %s because nJobs = %i", tmp_output_directory, options.nJobs )
    else:
        logger.info( "Output directory %s exists. Deleting.", tmp_output_directory )
        shutil.rmtree(tmp_output_directory)

try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(tmp_output_directory)
    logger.info( "Created output directory %s.", tmp_output_directory )
except:
    pass

# Anticipate & check output file
filename, ext = os.path.splitext( os.path.join(tmp_output_directory, sample.name + '.root') )
outfilename   = filename+ext

if not options.overwrite:
    if os.path.isfile(outfilename):
        logger.info( "Output file %s found.", outfilename)
        if checkRootFile( outfilename, checkForObjects=["Events"] ) and deepCheckRootFile( outfilename ) and deepCheckWeight( outfilename ):
            logger.info( "File already processed. Source: File check ok! Skipping." ) # Everything is fine, no overwriting
            sys.exit(0)
        else:
            logger.info( "File corrupt. Removing file from target." )
            os.remove( outfilename )
            logger.info( "Reprocessing." )
    else:
        logger.info( "Sample not processed yet." )
        logger.info( "Processing." )
else:
    logger.info( "Overwriting.")

# relocate original
sample.copy_files( os.path.join(tmp_output_directory, "input") )

treeFormulas = {}
if options.triggerSelection and isTriLep:
    # Trigger selection
    from tWZ.Tools.triggerSelector import triggerSelector
    ts           = triggerSelector(options.year)
    triggerCond  = ts.getSelection(options.samples[0] if isData else "MC", triggerList = ts.getTriggerList(sample) )
    treeFormulas["triggerDecision"] =  {'string':triggerCond}
    if isData:
        logger.info("Sample will have the following trigger skim: %s"%triggerCond)
        skimConds.append( triggerCond )


# turn on all branches to be flexible for filter cut in skimCond etc.
sample.chain.SetBranchStatus("*",1)

# this is the global selectionString
selectionString = '&&'.join(skimConds)

################################################################################
# top pt reweighting
from tWZ.Tools.topPtReweighting import getUnscaledTopPairPtReweightungFunction, getTopPtDrawString, getTopPtsForReweighting
# Decision based on sample name -> whether TTJets or TTLep is in the sample name
isTT = sample.name.startswith("TTJets") or sample.name.startswith("TTLep") or sample.name.startswith("TT_pow")
doTopPtReweighting = isTT and not options.noTopPtReweighting

if sample.name.startswith("TTLep"):
    sample.topScaleF = 1.002 ## found to be universal for years 2016-2018, and in principle negligible

if doTopPtReweighting:
    logger.info( "Sample will have top pt reweighting." )
    topPtReweightingFunc = getUnscaledTopPairPtReweightungFunction(selection = "dilep")
    # Compute x-sec scale factor on unweighted events
    if hasattr(sample, "topScaleF"):
        # If you don't want to get the SF for each subjob run the script and add the topScaleF to the sample
        topScaleF = sample.topScaleF
    else:
        reweighted  = sample.getYieldFromDraw( selectionString = selectionString, weightString = getTopPtDrawString(selection = "dilep") + '*genWeight')
        central     = sample.getYieldFromDraw( selectionString = selectionString, weightString = 'genWeight')

        topScaleF = central['val']/reweighted['val']

    logger.info( "Found topScaleF %f", topScaleF )
else:
    topScaleF = 1
    logger.info( "Sample will NOT have top pt reweighting. topScaleF=%f",topScaleF )

################################################################################
# CR reweighting 
if options.doCRReweighting:
    from tWZ.Tools.colorReconnectionReweighting import getCRWeight, getCRDrawString
    logger.info( "Sample will have CR reweighting." )
    #norm = sample.getYieldFromDraw( selectionString = selectionString, weightString = "genWeight" )
    norm = float(sample.chain.GetEntries(selectionString))
    CRScaleF = sample.getYieldFromDraw( selectionString = selectionString, weightString = getCRDrawString() )
    CRScaleF = CRScaleF['val']/norm#['val']
    logger.info(" Using norm: %s"%norm )
    logger.info(" Found CRScaleF: %s"%CRScaleF)
else:
    CRScaleF = 1
    logger.info( "Sample will NOT have CR reweighting. CRScaleF=%f",CRScaleF )

################################################################################
# branches to be kept for data and MC
branchKeepStrings_DATAMC = [\
    "run", "luminosityBlock", "event", "fixedGridRhoFastjetAll", "PV_npvs", "PV_npvsGood",
    "MET_*",
    "CaloMET_*",
    "RawMET_phi", "RawMET_pt", "RawMET_sumEt",
    "Flag_*",
    "nJet", "Jet_*",
    "nElectron", "Electron_*",
    "nMuon", "Muon_*",
    "L1PreFiringWeight_*",
]
branchKeepStrings_DATAMC += ["HLT_*", "PV_*"]

# NOT FOR UL?
if options.year == "2017":
    branchKeepStrings_DATAMC += [ "METFixEE2017_*" ]

#branches to be kept for MC samples only
branchKeepStrings_MC = [ "Generator_*", "GenPart_*", "nGenPart", "genWeight", "Pileup_nTrueInt","GenMET_*", "nGenJet", "GenJet_*", "genTtbarId"]
branchKeepStrings_MC.extend([ "*LHEScaleWeight", "*LHEPdfWeight", "LHEWeight_originalXWGTUP"])

#branches to be kept for data only
branchKeepStrings_DATA = [ ]

# Jet variables to be read from chain
jetCorrInfo = []
jetMCInfo   = ['genJetIdx/I','hadronFlavour/I', 'partonFlavour/I']

if isData:
    lumiScaleFactor=None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_DATA
    from FWCore.PythonUtilities.LumiList import LumiList
    # Apply golden JSON
    lumiList = LumiList(os.path.expandvars(sample.json))
    logger.info( "Loaded json %s", sample.json )
else:
    lumiScaleFactor = xSection*targetLumi/float(sample.normalization) if xSection is not None else None
    branchKeepStrings = branchKeepStrings_DATAMC + branchKeepStrings_MC

jetVars         = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F', 'pt_nom/F', 'corr_JER/F'] + jetCorrInfo
if isMC:
    jetVars     += jetMCInfo
    jetVars     += ['pt_jesTotalUp/F', 'pt_jesTotalDown/F', 'pt_jerUp/F', 'pt_jerDown/F', 'corr_JER/F', 'corr_JEC/F']
jetVarNames     = [x.split('/')[0] for x in jetVars]
genLepVars      = ['pt/F', 'phi/F', 'eta/F', 'pdgId/I', 'genPartIdxMother/I', 'status/I', 'statusFlags/I'] # some might have different types
genLepVarNames  = [x.split('/')[0] for x in genLepVars]
# those are for writing leptons
lepVars         = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F','pfRelIso03_all/F','mvaFall17V2Iso_WP90/O', 'mvaTTH/F', 'sip3d/F','lostHits/I','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','eleIndex/I','muIndex/I']
lepVarNames     = [x.split('/')[0] for x in lepVars]

read_variables = map(TreeVariable.fromString, [ 'MET_pt/F', 'MET_phi/F', 'run/I', 'luminosityBlock/I', 'event/l', 'PV_npvs/I', 'PV_npvsGood/I'] )
if options.year == "2017":
    read_variables += map(TreeVariable.fromString, [ 'METFixEE2017_pt/F', 'METFixEE2017_phi/F', 'METFixEE2017_pt_nom/F', 'METFixEE2017_phi_nom/F'])
    if isMC:
        read_variables += map(TreeVariable.fromString, [ 'METFixEE2017_pt_jesTotalUp/F', 'METFixEE2017_pt_jesTotalDown/F', 'METFixEE2017_pt_jerUp/F', 'METFixEE2017_pt_jerDown/F', 'METFixEE2017_pt_unclustEnDown/F', 'METFixEE2017_pt_unclustEnUp/F', 'METFixEE2017_phi_jesTotalUp/F', 'METFixEE2017_phi_jesTotalDown/F', 'METFixEE2017_phi_jerUp/F', 'METFixEE2017_phi_jerDown/F', 'METFixEE2017_phi_unclustEnDown/F', 'METFixEE2017_phi_unclustEnUp/F', 'METFixEE2017_pt_jer/F', 'METFixEE2017_phi_jer/F'])
# else:
    # These variables don't exist anymore in UL, MET variations are calculated with JMECorrector
    # read_variables += map(TreeVariable.fromString, [ 'MET_pt_nom/F', 'MET_phi_nom/F' ])
    # if isMC:
    #     read_variables += map(TreeVariable.fromString, [ 'MET_pt_jesTotalUp/F', 'MET_pt_jesTotalDown/F', 'MET_pt_jerUp/F', 'MET_pt_jerDown/F', 'MET_pt_unclustEnDown/F', 'MET_pt_unclustEnUp/F', 'MET_phi_jesTotalUp/F', 'MET_phi_jesTotalDown/F', 'MET_phi_jerUp/F', 'MET_phi_jerDown/F', 'MET_phi_unclustEnDown/F', 'MET_phi_unclustEnUp/F', 'MET_pt_jer/F', 'MET_phi_jer/F'])
if isMC:
    read_variables += map(TreeVariable.fromString, [ 'GenMET_pt/F', 'GenMET_phi/F' ])

new_variables = [ 'weight/F', 'year/I']

if options.triggerSelection and isTriLep:
    new_variables+= ['triggerDecision/I']

read_variables.append( TreeVariable.fromString('L1PreFiringWeight_Dn/F') )
read_variables.append( TreeVariable.fromString('L1PreFiringWeight_Nom/F') )
read_variables.append( TreeVariable.fromString('L1PreFiringWeight_Up/F') )

if isMC:
    read_variables += [ TreeVariable.fromString('Pileup_nTrueInt/F') ]
    # reading gen particles for top pt reweighting
    read_variables.append( TreeVariable.fromString('nGenPart/I') )
    read_variables.append( VectorTreeVariable.fromString('GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]', nMax=200 )) # default nMax is 100, which would lead to corrupt values in this case
    read_variables.append( TreeVariable.fromString('genWeight/F') )
    read_variables.append( TreeVariable.fromString('nGenJet/I') )
    read_variables.append( VectorTreeVariable.fromString('GenJet[pt/F,eta/F,phi/F,partonFlavour/I,hadronFlavour/I]' ) )
    new_variables.extend([ 'reweightTopPt/F', 'reweightPU/F','reweightPUUp/F','reweightPUDown/F', 'reweightL1Prefire/F', 'reweightL1PrefireUp/F', 'reweightL1PrefireDown/F'])
    if options.doCRReweighting:
        new_variables.append('reweightCR/F')

read_variables += [\
    TreeVariable.fromString('nElectron/I'),
    VectorTreeVariable.fromString('Electron[pt/F,eta/F,phi/F,pdgId/I,cutBased/I,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,lostHits/b,mvaFall17V2Iso_WP80/O,mvaFall17V2Iso_WP90/O,convVeto/O,dxy/F,dz/F,charge/I,deltaEtaSC/F,vidNestedWPBitmap/I,mvaTTH/F]'),
    TreeVariable.fromString('nMuon/I'),
    VectorTreeVariable.fromString('Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,dxy/F,dz/F,charge/I,mvaTTH/F]'), #jetNDauChargedMVASel/F
    TreeVariable.fromString('nJet/I'),
    VectorTreeVariable.fromString('Jet[%s]'% ( ','.join(jetVars) ) ) ]
if addReweights:
    read_variables.extend( ["nLHEReweightingWeight/I" ] )#, "nLHEReweighting/I", "LHEReweighting[Weight/F]" ] ) # need to set alias later

new_variables += [\
    'overlapRemoval/I','nlep/I',
    'JetGood[%s]'% ( ','.join(jetVars+['index/I']) + ',genPt/F' ),
    'met_pt/F', 'met_phi/F',
]

if sample.isData: new_variables.extend( ['jsonPassed/I','isData/I'] )
new_variables.extend( ['nBTag/I', 'm3/F', 'minDLmass/F'] )

new_variables.append( 'lep[%s]'% ( ','.join(lepVars) + ',mvaTOP/I' + ',mvaTOPWP/I' ) )
new_variables.append( 'Electron[mvaTOP/F]' )
new_variables.append( 'Muon[mvaTOP/F]' )

if isTriLep or isDiLep or isSingleLep:
    new_variables.extend( ['nGoodMuons/I', 'nGoodElectrons/I', 'nGoodLeptons/I' ] )
    new_variables.extend( ['l1_pt/F', 'l1_mvaTOP/F', 'l1_mvaTOPWP/I', 'l1_eta/F', 'l1_phi/F', 'l1_pdgId/I', 'l1_index/I', 'l1_jetPtRelv2/F', 'l1_jetPtRatiov2/F', 'l1_miniRelIso/F', 'l1_relIso03/F', 'l1_dxy/F', 'l1_dz/F', 'l1_mIsoWP/I', 'l1_eleIndex/I', 'l1_muIndex/I' ] )
    new_variables.extend( ['mlmZ_mass/F'])
    if isMC:
        new_variables.extend(['reweightLeptonSF/F', 'reweightLeptonSFUp/F', 'reweightLeptonSFDown/F'])
if isTriLep or isDiLep:
    new_variables.extend( ['l2_pt/F', 'l2_mvaTOP/F', 'l2_mvaTOPWP/I', 'l2_eta/F', 'l2_phi/F', 'l2_pdgId/I', 'l2_index/I', 'l2_jetPtRelv2/F', 'l2_jetPtRatiov2/F', 'l2_miniRelIso/F', 'l2_relIso03/F', 'l2_dxy/F', 'l2_dz/F', 'l2_mIsoWP/I', 'l2_eleIndex/I', 'l2_muIndex/I' ] )
    if isMC: new_variables.extend( \
        [   'genZ1_pt/F', 'genZ1_eta/F', 'genZ1_phi/F',
            'genZ2_pt/F', 'genZ1_eta/F', 'genZ1_phi/F',
            'reweightTrigger/F', 'reweightTriggerUp/F', 'reweightTriggerDown/F',
            'reweightLeptonTrackingSF/F',
         ] )
if isTriLep:
    new_variables.extend( ['l3_pt/F', 'l3_mvaTOP/F', 'l3_mvaTOPWP/I', 'l3_eta/F', 'l3_phi/F', 'l3_pdgId/I', 'l3_index/I', 'l3_jetPtRelv2/F', 'l3_jetPtRatiov2/F', 'l3_miniRelIso/F', 'l3_relIso03/F', 'l3_dxy/F', 'l3_dz/F', 'l3_mIsoWP/I', 'l3_eleIndex/I', 'l3_muIndex/I' ] )

if addReweights:
#    sample.chain.SetAlias("nLHEReweighting",       "nLHEReweightingWeight")
#    sample.chain.SetAlias("LHEReweighting_Weight", "LHEReweightingWeight")
    new_variables.append( TreeVariable.fromString("p[C/F]") )
    new_variables[-1].nMax = HyperPoly.get_ndof(weightInfo.nvar, options.interpolationOrder)
    new_variables.append( "chi2_ndof/F" )

## ttZ related variables
new_variables.extend( ['Z1_l1_index/I', 'Z1_l2_index/I', 'Z2_l1_index/I', 'Z2_l2_index/I', 'nonZ1_l1_index/I', 'nonZ1_l2_index/I'] )
for i in [1,2]:
    new_variables.extend( ['Z%i_pt/F'%i, 'Z%i_eta/F'%i, 'Z%i_phi/F'%i, 'Z%i_lldPhi/F'%i, 'Z%i_lldR/F'%i,  'Z%i_mass/F'%i, 'Z%i_cosThetaStar/F'%i] )

if options.checkTTGJetsOverlap:
    new_variables.extend( ['TTGJetsEventType/I'] )

if addSystematicVariations:
    for var in ['jesTotalUp', 'jesTotalDown', 'jerUp', 'jer', 'jerDown', 'unclustEnUp', 'unclustEnDown']:
        if not var.startswith('unclust'):
            new_variables.extend( ['nJetGood_'+var+'/I', 'nBTag_'+var+'/I'] )
        # new_variables.extend( ['met_pt_'+var+'/F', 'met_phi_'+var+'/F'] ) # MET variations are calculated with JMECorrector


# Btag weights Method 1a
for var in btagEff.btagWeightNames:
    if var!='MC':
        new_variables.append('reweightBTag_'+var+'/F')

if not options.skipNanoTools:
    ### nanoAOD postprocessor
    from importlib import import_module
    from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor   import PostProcessor
    from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel       import Collection
    from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop       import Module

    ## modules for nanoAOD postprocessor

    logger.info("Preparing nanoAOD postprocessing")
    logger.info("Will put files into directory %s", tmp_output_directory)

    if options.overwriteJEC is not None:
        JEC = options.overwriteJEC

    logger.info("Using JERs for MET significance")

    from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
    METBranchName = 'MET' if not options.year == "2017" else 'METFixEE2017'

    # check if files are available (e.g. if dpm is broken this should result in an error)
    for f in sample.files:
        if not checkRootFile(f):
            raise IOError ("File %s not available"%f)

    # remove empty files. this is necessary in 2018 because empty miniAOD files exist.
    sample.files = [ f for f in sample.files if nonEmptyFile(f) ]
    newFileList = []

    runPeriod = None
    if sample.isData:
        runString = sample.name.split('_')[1]
        assert str(options.year) in runString, "Could not obtain run period from sample name %s" % sample.name
        runPeriod = runString[-1]

    logger.info("Starting nanoAOD postprocessing")
    for f in sample.files:
        JMECorrector = createJMECorrector(
            isMC        = (not sample.isData),
            dataYear    = options.year,
            runPeriod   = runPeriod,
            jesUncert   = "Total",
            jetType     = "AK4PFchs",
            metBranchName = METBranchName,
            isFastSim   = False,
            applySmearing = False)

        modules = [ JMECorrector() ]

        # need a hash to avoid data loss
        file_hash = str(hash(f))
        p = PostProcessor(tmp_output_directory, [f], cut=selectionString, modules=modules, postfix="_for_%s_%s"%(sample.name, file_hash))
        if not options.reuseNanoAOD:
            p.run()
        newFileList += [tmp_output_directory + '/' + f.split('/')[-1].replace('.root', '_for_%s_%s.root'%(sample.name, file_hash))]
    logger.info("Done. Replacing input files for further processing.")

    sample.files = newFileList
    sample.clear()

# Define a reader

logger.info( "Running with selectionString %s", selectionString )
reader = sample.treeReader( \
    variables = read_variables,
    selectionString = selectionString
    )

# using miniRelIso 0.2 as baseline
#eleSelector_ = eleSelector( "CBtight", year = options.year )
#muSelector_  = muonSelector("medium",  year = options.year )
eleSelector_ = eleSelector( "mvaTOPVL", year = options.year )
muSelector_  = muonSelector("mvaTOPVL", year = options.year )

################################################################################
# FILL EVENT INFO HERE
################################################################################
def filler( event ):
    # shortcut
    r = reader.event
    #workaround  = (r.run, r.luminosityBlock, r.event) # some fastsim files seem to have issues, apparently solved by this.
    event.isData = s.isData
    event.year   = yearint
    event.overlapRemoval = 1

    if isMC:

        # GEN Particles
        gPart = getGenPartsAll(r)
        # GEN Jets
        gJets = getJets( r, jetVars=['pt','eta','phi','mass','partonFlavour','hadronFlavour','index'], jetColl="GenJet" )

        genLepsFromZ    = genLepFromZ(gPart)
        genZs           = getSortedZCandidates(genLepsFromZ)
        if len(genZs)>0:
            event.genZ_mass = genZs[0][0]

    if isMC:
        if hasattr(r, "genWeight"):
            event.weight = lumiScaleFactor*r.genWeight if lumiScaleFactor is not None else 1
        else:
            event.weight = lumiScaleFactor if lumiScaleFactor is not None else 1
    elif sample.isData:
        event.weight = 1
    else:
        raise NotImplementedError( "isMC %r isData %r " % (isMC, isData) )

    # lumi lists and vetos
    if sample.isData:
        #event.vetoPassed  = vetoList.passesVeto(r.run, r.lumi, r.evt)
        event.jsonPassed  = lumiList.contains(r.run, r.luminosityBlock)
        # apply JSON to data via event weight
        if not event.jsonPassed: event.weight=0
        # store decision to use after filler has been executed
        event.jsonPassed_ = event.jsonPassed

    ################################################################################
    # PU reweighting
    if isMC and hasattr(r, "Pileup_nTrueInt"):
        from Analysis.Tools.puWeightsUL			import getPUReweight
        event.reweightPU     = getPUReweight( r.Pileup_nTrueInt, year=options.year, weight="nominal") 
        event.reweightPUDown = getPUReweight( r.Pileup_nTrueInt, year=options.year, weight="down" )
        event.reweightPUUp   = getPUReweight( r.Pileup_nTrueInt, year=options.year, weight="up" )

    ################################################################################
    # top pt reweighting
    if isMC:
        event.reweightTopPt     = topPtReweightingFunc(getTopPtsForReweighting(r)) * topScaleF if doTopPtReweighting else 1.

    ################################################################################
    # Trigger Decision
    if options.triggerSelection and isTriLep:
        event.triggerDecision = int(treeFormulas['triggerDecision']['TTreeFormula'].EvalInstance())

    if addReweights:

        include_missing_refpoint = False
        if weightInfo.nid == r.nLHEReweightingWeight + 1:
            logger.debug( "addReweights: pkl file has %i base-points but nLHEReweightWeight=%i, hence it is likely that the ref-point is among the base points and is missing. Fingers crossed.", weightInfo.nid, r.nLHEReweightingWeight )
            if ref_point_index is None:
                raise RuntimeError( "weightInfo.nid == r.nLHEReweightingWeight + 1 but ref_point_index = None -> something is wrong." )
            include_missing_refpoint = True
        elif weightInfo.nid == r.nLHEReweightingWeight:
            pass
        else:
            raise RuntimeError("reweight_pkl and nLHEReweightWeight are inconsistent.")

        # here we check the consistency with miniAOD, hence multiply with LHEWeight_originalXWGTUP
        weights = [reader.sample.chain.GetLeaf("LHEReweightingWeight").GetValue(i_weight) for i_weight in range(r.nLHEReweightingWeight)]
        if include_missing_refpoint:
            weights = weights[:ref_point_index] + [1] + weights[ref_point_index:]

        coeff           = hyperPoly.get_parametrization( weights )
        event.np        = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof( coeff, weights )
        #print event.chi2_ndof, event.np, weights, coeff
        if event.chi2_ndof > 10**-6:
            logger.warning( "chi2_ndof is large: %f", event.chi2_ndof )

        for n in xrange( hyperPoly.ndof ):
            event.p_C[n] = coeff[n]
    #    p_C_nanoAOD.append( [ coeff[n] for n in xrange(hyperPoly.ndof) ] )

    allSlimmedJets      = getJets(r)
    event.reweightL1Prefire, event.reweightL1PrefireUp, event.reweightL1PrefireDown = r.L1PreFiringWeight_Nom, r.L1PreFiringWeight_Up, r.L1PreFiringWeight_Dn

    # get leptons before jets in order to clean jets
    electrons_pt10  = getGoodElectrons(r, ele_selector = eleSelector_)
    muons_pt10      = getGoodMuons    (r, mu_selector  = muSelector_ )

    for e in electrons_pt10:
        e['pdgId']      = int( -11*e['charge'] )
        e['eleIndex']   = e['index']
        e['muIndex']    = -1
        e['mvaTOP']     = e['mvaTTH']
    for m in muons_pt10:
        m['pdgId']      = int( -13*m['charge'] )
        m['muIndex']    = m['index']
        m['eleIndex']   = -1
        m['mvaTOP']     = m['mvaTTH']
        # print m['jetNDauChargedMVASel']

    # make list of leptons
    leptons = electrons_pt10+muons_pt10
    leptons.sort(key = lambda p:-p['pt'])

    for iLep, lep in enumerate(leptons):
        lep['index']    = iLep
        lep['mvaTOP']   = lep['mvaTTH']
        lep['mvaTOPWP'] = mvaTopWP(lep['mvaTOP'], lep['pdgId'])

    fill_vector_collection( event, "lep", lepVarNames, leptons)
    event.nlep = len(leptons)

    event.minDLmass = getMinDLMass(leptons)

    # now get jets, cleaned against good leptons
    all_jets     = getJets(r, jetVars=jetVarNames)
    clean_jets,_ = cleanJetsAndLeptons( all_jets, leptons )
    clean_jets_acc = filter(lambda j:abs(j['eta'])<2.4, clean_jets)
    jets         = filter(lambda j:j['pt']>30, clean_jets_acc)
    bJets        = filter(lambda j:isBJet(j, tagger=b_tagger, year=options.year) and abs(j['eta'])<=2.4    , jets)
    nonBJets     = filter(lambda j:not ( isBJet(j, tagger=b_tagger, year=options.year) and abs(j['eta'])<=2.4 ), jets)

    # store the correct MET (EE Fix for 2017)

    if options.year == 2017:# and not options.fastSim:
        # v2 recipe. Could also use our own recipe
        event.met_pt    = r.METFixEE2017_pt_nom
        event.met_phi   = r.METFixEE2017_phi_nom
    else:
        event.met_pt    = r.MET_pt
        event.met_phi   = r.MET_phi

    # Filling jets
    maxNJet = 100
    store_jets = jets #if not options.keepAllJets else soft_jets + jets
    store_jets = store_jets[:maxNJet]
    store_jets.sort( key = lambda j:-j['pt'])
    event.nJetGood   = len(store_jets)
    for iJet, jet in enumerate(store_jets):
        event.JetGood_index[iJet] = jet['index']
        for b in jetVarNames:
            getattr(event, "JetGood_"+b)[iJet] = jet[b]
        if isMC:
            if store_jets[iJet]['genJetIdx'] >= 0:
                if r.nGenJet<maxNJet:
                    try:
                        event.JetGood_genPt[iJet] = r.GenJet_pt[store_jets[iJet]['genJetIdx']]
                    except IndexError:
                        event.JetGood_genPt[iJet] = -1
                else:
                    event.JetGood_genPt[iJet] = -1
        getattr(event, "JetGood_pt")[iJet] = jet['pt']

    if isMC and options.doCRReweighting:
        event.reweightCR = getCRWeight(event.nJetGood)

    event.nBTag      = len(bJets)
    event.m3, _,_,_  = m3(jets)

    # we got all objects, so let's make dR

    jets_sys      = {}
    bjets_sys     = {}
    nonBjets_sys  = {}

    if addSystematicVariations:
        for var in ['jesTotalUp', 'jesTotalDown', 'jerUp', 'jerDown', 'unclustEnUp', 'unclustEnDown']: # don't use 'jer' as of now
            # MET variations are calculated with JMECorrector, not here
            # setattr(event, 'met_pt_'+var,  getattr(r, 'METFixEE2017_pt_'+var)  if options.year == 2017 else getattr(r, 'MET_pt_'+var) )
            # setattr(event, 'met_phi_'+var, getattr(r, 'METFixEE2017_phi_'+var) if options.year == 2017 else getattr(r, 'MET_phi_'+var) )
            if not var.startswith('unclust'):
                corrFactor = 'corr_JER' if var == 'jer' else None
                jets_sys[var]       = filter(lambda j:j['pt_'+var]>30, clean_jets_acc)
                bjets_sys[var]      = filter(lambda j: isBJet(j) and abs(j['eta'])<2.4, jets_sys[var])
                nonBjets_sys[var]   = filter(lambda j: not ( isBJet(j) and abs(j['eta'])<2.4), jets_sys[var])

                setattr(event, "nJetGood_"+var, len(jets_sys[var]))
                setattr(event, "nBTag_"+var,    len(bjets_sys[var]))

    if isSingleLep or isTriLep or isDiLep:
        event.nGoodMuons      = len(filter( lambda l:abs(l['pdgId'])==13, leptons))
        event.nGoodElectrons  = len(filter( lambda l:abs(l['pdgId'])==11, leptons))
        event.nGoodLeptons    = len(leptons)

        if len(leptons)>=1:
            event.l1_pt         = leptons[0]['pt']
            event.l1_mvaTOP     = leptons[0]['mvaTOP']
            event.l1_mvaTOPWP   = leptons[0]['mvaTOPWP']
            event.l1_eta        = leptons[0]['eta']
            event.l1_phi        = leptons[0]['phi']
            event.l1_pdgId      = leptons[0]['pdgId']
            event.l1_index      = leptons[0]['index']
            event.l1_miniRelIso = leptons[0]['miniPFRelIso_all']
            event.l1_relIso03   = leptons[0]['pfRelIso03_all']
            event.l1_dxy        = leptons[0]['dxy']
            event.l1_dz         = leptons[0]['dz']
            event.l1_eleIndex   = leptons[0]['eleIndex']
            event.l1_muIndex    = leptons[0]['muIndex']

        # For TTZ studies: find Z boson candidate, and use third lepton to calculate mt
        (event.mlmZ_mass, zl1, zl2) = closestOSDLMassToMZ(leptons)

        if isMC:
            trig_eff, trig_eff_err =  triggerEff.getSF(leptons)
            event.reweightTrigger       = trig_eff
            event.reweightTriggerUp     = trig_eff + trig_eff_err
            event.reweightTriggerDown   = trig_eff - trig_eff_err

            leptonsForSF   = ( leptons[:2] if isDiLep else (leptons[:3] if isTriLep else leptons[:1]) )
            #leptonSFValues = [ leptonSF.getSF(pdgId=l['pdgId'], pt=l['pt'], eta=((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta'])) for l in leptonsForSF ]
            #event.reweightLeptonSF     = reduce(mul, [sf[0] for sf in leptonSFValues], 1)
            #event.reweightLeptonSFDown = reduce(mul, [sf[1] for sf in leptonSFValues], 1)
            #event.reweightLeptonSFUp   = reduce(mul, [sf[2] for sf in leptonSFValues], 1)
            #if event.reweightLeptonSF ==0:
            #    logger.error( "reweightLeptonSF is zero!")

            #event.reweightLeptonTrackingSF   = reduce(mul, [leptonTrackingSF.getSF(pdgId = l['pdgId'], pt = l['pt'], eta = ((l['eta'] + l['deltaEtaSC']) if abs(l['pdgId'])==11 else l['eta']))  for l in leptonsForSF], 1)

    if isTriLep or isDiLep:
        if len(leptons)>=2:
            event.l2_pt         = leptons[1]['pt']
            event.l2_mvaTOP     = leptons[1]['mvaTOP']
            event.l2_mvaTOPWP   = leptons[1]['mvaTOPWP']
            event.l2_eta        = leptons[1]['eta']
            event.l2_phi        = leptons[1]['phi']
            event.l2_pdgId      = leptons[1]['pdgId']
            event.l2_index      = leptons[1]['index']
            event.l2_miniRelIso = leptons[1]['miniPFRelIso_all']
            event.l2_relIso03   = leptons[1]['pfRelIso03_all']
            event.l2_dxy        = leptons[1]['dxy']
            event.l2_dz         = leptons[1]['dz']
            event.l2_eleIndex   = leptons[1]['eleIndex']
            event.l2_muIndex    = leptons[1]['muIndex']

            if isMC:
              genZs = getGenZs(gPart)
              if len(genZs)>0:
                  event.genZ1_pt  = genZs[0]['pt']
                  event.genZ1_eta = genZs[0]['eta']
                  event.genZ1_phi = genZs[0]['phi']
              else:
                  event.genZ1_pt  =float('nan')
                  event.genZ1_eta =float('nan')
                  event.genZ1_phi =float('nan')
              if len(genZs)>1:
                  event.genZ2_pt  = genZs[1]['pt']
                  event.genZ2_eta = genZs[1]['eta']
                  event.genZ2_phi = genZs[1]['phi']
              else:
                  event.genZ2_pt  =float('nan')
                  event.genZ2_eta =float('nan')
                  event.genZ2_phi =float('nan')

            # for quadlep stuff
            allZCands = getSortedZCandidates(leptons)
            Z_vectors = []
            for i in [0,1]:
                if len(allZCands) > i:
                    (Z_mass, Z_l1_tightLepton_index, Z_l2_tightLepton_index) = allZCands[i]
                    Z_l1_index = Z_l1_tightLepton_index if Z_l1_tightLepton_index>=0 else -1
                    Z_l2_index = Z_l2_tightLepton_index if Z_l2_tightLepton_index>=0 else -1
                    setattr(event, "Z%i_mass"%(i+1),       Z_mass)
                    setattr(event, "Z%i_l1_index"%(i+1),   Z_l1_index)
                    setattr(event, "Z%i_l2_index"%(i+1),   Z_l2_index)
                    Z_l1 = ROOT.TLorentzVector()
                    Z_l1.SetPtEtaPhiM(leptons[Z_l1_index]['pt'], leptons[Z_l1_index]['eta'], leptons[Z_l1_index]['phi'], 0 )
                    Z_l2 = ROOT.TLorentzVector()
                    Z_l2.SetPtEtaPhiM(leptons[Z_l2_index]['pt'], leptons[Z_l2_index]['eta'], leptons[Z_l2_index]['phi'], 0 )
                    Z = Z_l1 + Z_l2
                    setattr(event, "Z%i_pt"%(i+1),         Z.Pt())
                    setattr(event, "Z%i_eta"%(i+1),        Z.Eta())
                    setattr(event, "Z%i_phi"%(i+1),        Z.Phi())
                    setattr(event, "Z%i_lldPhi"%(i+1),     deltaPhi(Z_l1.Phi(), Z_l2.Phi()))
                    setattr(event, "Z%i_lldR"%(i+1),       deltaR(leptons[Z_l1_index], leptons[Z_l2_index]) )
                    if Z_l1_index>=0:
                        lm_Z_index = Z_l1_index if event.lep_pdgId[Z_l1_index] > 0 else Z_l2_index
                        setattr(event, "Z%i_cosThetaStar"%(i+1), cosThetaStar(Z_mass, Z.Pt(), Z.Eta(), Z.Phi(), event.lep_pt[lm_Z_index], event.lep_eta[lm_Z_index], event.lep_phi[lm_Z_index]) )

                    Z_vectors.append(Z)

            # take the leptons that are not from the leading Z candidate and assign them as nonZ, ignorant about if they actually from a Z candidate
            # As a start, take the leading two leptons as non-Z. To be overwritten as soon as we have a Z candidate, otherwise one lepton can be both from Z and non-Z
            if len(leptons)>0:
                event.nonZ1_l1_index = leptons[0]['index']
            if len(leptons)>1:
                event.nonZ1_l2_index = leptons[1]['index']
            if len(allZCands)>0:
                # reset nonZ1_leptons
                event.nonZ1_l1_index = -1
                event.nonZ1_l2_index = -1
                nonZ1_tightLepton_indices = [ i for i in range(len(leptons)) if i not in [allZCands[0][1], allZCands[0][2]] ]

                event.nonZ1_l1_index = nonZ1_tightLepton_indices[0] if len(nonZ1_tightLepton_indices)>0 else -1
                event.nonZ1_l2_index = nonZ1_tightLepton_indices[1] if len(nonZ1_tightLepton_indices)>1 else -1

        if len(leptons)>=3:
            event.l3_pt         = leptons[2]['pt']
            event.l3_mvaTOP     = leptons[2]['mvaTOP']
            event.l3_mvaTOPWP   = leptons[2]['mvaTOPWP']
            event.l3_eta        = leptons[2]['eta']
            event.l3_phi        = leptons[2]['phi']
            event.l3_pdgId      = leptons[2]['pdgId']
            event.l3_index      = leptons[2]['index']
            event.l3_miniRelIso = leptons[2]['miniPFRelIso_all']
            event.l3_relIso03   = leptons[2]['pfRelIso03_all']
            event.l3_dxy        = leptons[2]['dxy']
            event.l3_dz         = leptons[2]['dz']
            event.l3_eleIndex   = leptons[2]['eleIndex']
            event.l3_muIndex    = leptons[2]['muIndex']

    #if addSystematicVariations:
    # B tagging weights method 1a
    if isMC:
        for j in jets:
            btagEff.addBTagEffToJet(j)
        for var in btagEff.btagWeightNames:
            if var!='MC':
                setattr(event, 'reweightBTag_'+var, btagEff.getBTagSF_1a( var, bJets, filter( lambda j: abs(j['eta'])<2.4, nonBJets ) ) )


# Create a maker. Maker class will be compiled. This instance will be used as a parent in the loop
treeMaker_parent = TreeMaker(
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) if type(x)==type("") else x for x in new_variables ],
    treeName = "Events"
    )

# Split input in ranges
eventRanges = reader.getEventRanges( maxNEvents = options.eventsPerJob, minJobs = options.minNJobs )

logger.info( "Splitting into %i ranges of %i events on average. FileBasedSplitting: %s. Job number %s",
        len(eventRanges),
        (eventRanges[-1][1] - eventRanges[0][0])/len(eventRanges),
        'Yes',
        options.job)

#Define all jobs
jobs = [(i, eventRanges[i]) for i in range(len(eventRanges))]

filename, ext = os.path.splitext( os.path.join(tmp_output_directory, sample.name + '.root') )

if len(eventRanges)>1:
    raise RuntimeError("Using fileBasedSplitting but have more than one event range!")

clonedEvents = 0
convertedEvents = 0
outputLumiList = {}
for ievtRange, eventRange in enumerate( eventRanges ):

    logger.info( "Processing range %i/%i from %i to %i which are %i events.",  ievtRange, len(eventRanges), eventRange[0], eventRange[1], eventRange[1]-eventRange[0] )

    _logger.   add_fileHandler( outfilename.replace('.root', '.log'), options.logLevel )
    _logger_rt.add_fileHandler( outfilename.replace('.root', '_rt.log'), options.logLevel )

    tmp_gdirectory = ROOT.gDirectory
    outputfile = ROOT.TFile.Open(outfilename, 'recreate')
    tmp_gdirectory.cd()

    if options.small:
        logger.info("Running 'small'. Not more than 10000 events")
        nMaxEvents = eventRange[1]-eventRange[0]
        eventRange = ( eventRange[0], eventRange[0] +  min( [nMaxEvents, 10000] ) )

    # Set the reader to the event range
    reader.setEventRange( eventRange )

    clonedTree = reader.cloneTree( branchKeepStrings, newTreename = "Events", rootfile = outputfile )
    clonedEvents += clonedTree.GetEntries()

    # Add the TTreeFormulas
    for formula in treeFormulas.keys():
        treeFormulas[formula]['TTreeFormula'] = ROOT.TTreeFormula(formula, treeFormulas[formula]['string'], clonedTree )

    # Clone the empty maker in order to avoid recompilation at every loop iteration
    maker = treeMaker_parent.cloneWithoutCompile( externalTree = clonedTree )

    maker.start()
    # Do the thing
    reader.start()
    reader.sample.chain.SetBranchStatus("*",1)
    while reader.run():

        maker.run()
        if sample.isData:
            if maker.event.jsonPassed_:
                if reader.event.run not in outputLumiList.keys():
                    outputLumiList[reader.event.run] = set([reader.event.luminosityBlock])
                else:
                    if reader.event.luminosityBlock not in outputLumiList[reader.event.run]:
                        outputLumiList[reader.event.run].add(reader.event.luminosityBlock)

    convertedEvents += maker.tree.GetEntries()
    maker.tree.Write()
    outputfile.Close()
    logger.info( "Written %s", outfilename)

  # Destroy the TTree
    maker.clear()
    sample.clear()


logger.info( "Converted %i events of %i, cloned %i",  convertedEvents, reader.nEvents , clonedEvents )

# Storing JSON file of processed events
if sample.isData and convertedEvents>0: # avoid json to be overwritten in cases where a root file was found already
    jsonFile = filename+'_%s.json'%(0 if options.nJobs==1 else options.job)
    LumiList( runsAndLumis = outputLumiList ).writeJSON(jsonFile)
    logger.info( "Written JSON file %s", jsonFile )

if not ( options.keepNanoAOD or options.reuseNanoAOD) and not options.skipNanoTools:
    for f in sample.files:
        try:
            os.remove(f)
            logger.info("Removed nanoAOD file: %s", f)
        except OSError:
            logger.info("nanoAOD file %s seems to be not there", f)

logger.info("Copying log file to %s", storage_directory )
copyLog = subprocess.call(['cp', logFile, storage_directory] )
if copyLog:
    logger.info( "Copying log from %s to %s failed", logFile, storage_directory)
else:
    logger.info( "Successfully copied log file" )
    os.remove(logFile)
    logger.info( "Removed temporary log file" )

if checkRootFile( outfilename, checkForObjects=["Events"] ) and deepCheckRootFile( outfilename ) and deepCheckWeight( outfilename ):
    logger.info( "Target: File check ok!" )
else:
    logger.info( "Corrupt rootfile! Removing file: %s"%outfilename )
    os.remove( outfilename )

for item in os.listdir(tmp_output_directory):
    s = os.path.join(tmp_output_directory, item)
    if not os.path.isdir(s):
        shutil.copy(s, storage_directory)
logger.info( "Done copying to storage directory %s", storage_directory)

# close all log files before deleting the tmp directory
for logger_ in [logger, logger_rt]:
    for handler in logger_.handlers:
        handler.close()
        logger_.removeHandler(handler)

if os.path.exists(tmp_output_directory):
    shutil.rmtree(tmp_output_directory)
    logger.info( "Cleaned tmp directory %s", tmp_output_directory )

# There is a double free corruption due to stupid ROOT memory management which leads to a non-zero exit code
# Thus the job is resubmitted on condor even if the output is ok
# Current idea is that the problem is with xrootd having a non-closed root file
sample.clear()

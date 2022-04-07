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
from tWZ.Tools.helpers             import deltaR, deltaPhi, bestDRMatchInCollection, nonEmptyFile, cosThetaStar, m3
from tWZ.Tools.objectSelection     import getMuons, getElectrons, muonSelector, eleSelector, getGoodMuons, getGoodElectrons, isBJet, getGenPartsAll, getJets, mvaTopWP, getPhotons

from tWZ.Tools.leptonSF            import leptonSF as leptonSF_
from tWZ.Tools.mcTools import pdgToName, GenSearch, B_mesons, D_mesons, B_mesons_abs, D_mesons_abs
genSearch = GenSearch()

from Analysis.Tools.metFilters               import getFilterCut
# from Analysis.Tools.overlapRemovalTTG        import hasMesonMother, getParentIds
from Analysis.Tools.puProfileDirDB           import puProfile
from Analysis.Tools.L1PrefireWeight          import L1PrefireWeight
from Analysis.Tools.LeptonTrackingEfficiency import LeptonTrackingEfficiency
from Analysis.Tools.isrWeight                import ISRweight
from Analysis.Tools.helpers                  import checkRootFile, deepCheckRootFile, deepCheckWeight
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
    argParser.add_argument('--nJobs',       action='store',         nargs='?',  type=int, default=1,                                    help="Maximum number of simultaneous jobs." )
    argParser.add_argument('--job',         action='store',                     type=int, default=0,                                    help="Run only jobs i" )
    argParser.add_argument('--minNJobs',    action='store',         nargs='?',  type=int, default=1,                                    help="Minimum number of simultaneous jobs." )
    argParser.add_argument('--targetDir',   action='store',         nargs='?',  type=str, default=user.postprocessing_output_directory, help="Name of the directory the post-processed files will be saved" )
    argParser.add_argument('--processingEra', action='store',       nargs='?',  type=str, default='postProcessed_80X_v22',              help="Name of the processing era" )
    argParser.add_argument('--skim',        action='store',         nargs='?',  type=str, default='1jmu',                               help="Skim conditions to be applied for post-processing" )
    argParser.add_argument('--trigger',     action='store',         nargs='*',  type=str, default=['HLT_Mu3_PFJet40'],                  help="List of triggers" )
    argParser.add_argument('--LHEHTCut',    action='store',         nargs='?',  type=int, default=-1,                                   help="LHE cut." )
    argParser.add_argument('--year',        action='store',                     type=int,                                               help="Which year?" )
    argParser.add_argument('--overwriteJEC',action='store',                               default=None,                                 help="Overwrite JEC?" )
    argParser.add_argument('--overwrite',   action='store_true',                                                                        help="Overwrite existing output files, bool flag set to True  if used" )
    argParser.add_argument('--small',       action='store_true',                                                                        help="Run the file on a small sample (for test purpose), bool flag set to True if used" )
    argParser.add_argument('--flagTTGamma', action='store_true',                                                                        help="Is ttgamma?" )
    argParser.add_argument('--flagTTBar',   action='store_true',                                                                        help="Is ttbar?" )
    argParser.add_argument('--doCRReweighting',             action='store_true',                                                        help="color reconnection reweighting?")
    #argParser.add_argument('--skipGenLepMatching',          action='store_true',                                                        help="skip matched genleps??" )
    argParser.add_argument('--checkTTGJetsOverlap',         action='store_true',                                                        help="Keep TTGJetsEventType which can be used to clean TTG events from TTJets samples" )
    argParser.add_argument('--forceProxy',                  action='store_true',                                                        help="Don't check certificate")
    argParser.add_argument('--skipNanoTools',               action='store_true',                                                        help="Skipt the nanoAOD tools step for computing JEC/JER/MET etc uncertainties")
    argParser.add_argument('--keepNanoAOD',                 action='store_true',                                                        help="Keep nanoAOD output?")
    argParser.add_argument('--reuseNanoAOD',                action='store_true',                                                        help="Keep nanoAOD output?")
    argParser.add_argument('--reapplyJECS',                 action='store_true',                                                        help="Reapply JECs to data?")
    argParser.add_argument('--reduceSizeBy',                action='store',     type=int,                                               help="Reduce the size of the sample by a factor of...")
    argParser.add_argument('--event',                       action='store',     type=int, default=-1,                                   help="Just process event no")

    return argParser

options = get_parser().parse_args()

# for t in options.trigger:
#     print "{trigger}/O".format(trigger =t)

# print TreeVariable.fromString("HLT_Ele8_CaloIdM_TrackIdM_PFJet30/O")
# print "+"*90
# print (
#     [
#         TreeVariable.fromString("{trigger}/O".format(trigger =t)) for t in options.trigger
#     ]
# )
# exit(0)

# print ["TreeVariable.fromString('{trigger}/I')".format(trigger =t) for t in options.trigger]
# exit(0)

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
                #if type(obj[var]) == type("string"):
                #    obj[var] = int(ord(obj[var]))
                if type(obj[var]) == type(True):
                    obj[var] = int(obj[var])
                #print collection_name+"_"+var, getattr(event, collection_name+"_"+var), i_obj, var, obj[var]
                getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

# Skim condition
if options.skim.startswith('1j1mu'):
    skimConds = ["Sum$(Jet_pt>40)>=1&&Sum$(Muon_pt>3.5)>=1"]
    #triggerCond  = "HLT_Mu3_PFJet40" 
elif options.skim.startswith('1j1ele'):
    skimConds = ["Sum$(Jet_pt>40)>=1&&Sum$(Electron_pt>5)>=1"]
    #triggerCond  = "HLT_Ele8_CaloIdM_TrackIdM_PFJet30"
    #triggerCond  = "HLT_Ele8_CaloIdM_TrackIdM_PFJet30", "_HLT_Ele12_CaloIdM_TrackIdM_PFJet30", "_HLT_Ele17_CaloIdM_TrackIdM_PFJet30", "_HLT_Ele23_CaloIdM_TrackIdM_PFJet30" 

# trigger condition
triggerCond = '('+'||'.join(options.trigger)+')'
# root -l /eos/vbc/incoming///store/group/phys_susy/stops2l/topNanoAOD/v6-1-2/2016/JetHT/TopNanoAODv6-1-2-6_JetHT_Run2016B_ver2/200909_064328/0000/tree_857.root

#Samples: Load samples
maxNfiles = 1 if options.small else None
if options.small:
    options.job = 0
    #options.nJobs = 10000 # set high to just run over 1 input file

if options.year == 2016:
    from Samples.nanoAOD.Summer16_private_nanoAODv6         import allSamples as mcSamples
    from Samples.nanoAOD.Run2016_private_nanoAODv6          import allSamples as dataSamples
    allSamples = mcSamples + dataSamples
elif options.year == 2017:
    from Samples.nanoAOD.Fall17_private_nanoAODv6           import allSamples as mcSamples
    from Samples.nanoAOD.Run2017_private_nanoAODv6          import allSamples as dataSamples
    allSamples = mcSamples + dataSamples
elif options.year == 2018:
    from Samples.nanoAOD.Autumn18_private_nanoAODv6         import allSamples as mcSamples
    from Samples.nanoAOD.Run2018_private_nanoAODv6          import allSamples as dataSamples
    allSamples = mcSamples + dataSamples

samples = []
for selectedSamples in options.samples:
    for sample in allSamples:
        if selectedSamples == sample.name:
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

L1PW = L1PrefireWeight(options.year)

treeFormulas = {"triggerDecision": {'string':triggerCond} }
logger.info("Sample will have the following trigger skim: %s"%triggerCond)
skimConds.append( triggerCond )

# apply MET filter
skimConds.append( getFilterCut(options.year, isData=isData, ignoreJSON=True, skipWeight=True) )

isr                 = ISRweight()

#Samples: combine if more than one
if len(samples)>1:
    sample_name =  samples[0].name+"_comb"
    logger.info( "Combining samples %s to %s.", ",".join(s.name for s in samples), sample_name )
    sample      = Sample.combine(sample_name, samples, maxN = maxNfiles)
    sampleForPU = Sample.combine(sample_name, samples, maxN = -1)
elif len(samples)==1:
    sample      = samples[0]
    sampleForPU = samples[0]
    if options.small:
        sample.reduceFiles( to = maxNfiles )
else:
    raise ValueError( "Need at least one sample. Got %r",samples )


# print sample.files
# exit(0)

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

if isMC:
    from Analysis.Tools.puReweighting import getReweightingFunction
    if options.year == 2016:
        nTrueInt_puRW       = getReweightingFunction(data="PU_2016_35920_XSecCentral", mc="Summer16")
        nTrueInt_puRWDown   = getReweightingFunction(data="PU_2016_35920_XSecDown",    mc="Summer16")
        nTrueInt_puRWVDown  = getReweightingFunction(data="PU_2016_35920_XSecVDown",   mc="Summer16")
        nTrueInt_puRWUp     = getReweightingFunction(data="PU_2016_35920_XSecUp",      mc="Summer16")
        nTrueInt_puRWVUp    = getReweightingFunction(data="PU_2016_35920_XSecVUp",     mc="Summer16")
        nTrueInt_puRWVVUp   = getReweightingFunction(data="PU_2016_35920_XSecVVUp",    mc="Summer16")
    elif options.year == 2017:
        # keep the weight name for now. Should we update to a more general one?
        puProfiles = puProfile( source_sample = sampleForPU )
        mcHist = puProfiles.cachedTemplate( selection="( 1 )", weight='genWeight', overwrite=False ) # use genWeight for amc@NLO samples. No problems encountered so far
        nTrueInt_puRW       = getReweightingFunction(data="PU_2017_41530_XSecCentral",  mc=mcHist)
        nTrueInt_puRWDown   = getReweightingFunction(data="PU_2017_41530_XSecDown",     mc=mcHist)
        nTrueInt_puRWVDown  = getReweightingFunction(data="PU_2017_41530_XSecVDown",    mc=mcHist)
        nTrueInt_puRWUp     = getReweightingFunction(data="PU_2017_41530_XSecUp",       mc=mcHist)
        nTrueInt_puRWVUp    = getReweightingFunction(data="PU_2017_41530_XSecVUp",      mc=mcHist)
        nTrueInt_puRWVVUp   = getReweightingFunction(data="PU_2017_41530_XSecVVUp",     mc=mcHist)
    elif options.year == 2018:
        # keep the weight name for now. Should we update to a more general one?
        nTrueInt_puRW       = getReweightingFunction(data="PU_2018_59740_XSecCentral",  mc="Autumn18")
        nTrueInt_puRWDown   = getReweightingFunction(data="PU_2018_59740_XSecDown",     mc="Autumn18")
        nTrueInt_puRWVDown  = getReweightingFunction(data="PU_2018_59740_XSecVDown",    mc="Autumn18")
        nTrueInt_puRWUp     = getReweightingFunction(data="PU_2018_59740_XSecUp",       mc="Autumn18")
        nTrueInt_puRWVUp    = getReweightingFunction(data="PU_2018_59740_XSecVUp",      mc="Autumn18")
        nTrueInt_puRWVVUp   = getReweightingFunction(data="PU_2018_59740_XSecVVUp",     mc="Autumn18")

## lepton SFs
#leptonTrackingSF    = LeptonTrackingEfficiency(options.year)
#leptonSF            = leptonSF_(options.year)

options.skim = options.skim + '_small' if options.small else options.skim

# LHE cut (DY samples)
if options.LHEHTCut>0:
    sample.name+="_lheHT"+str(options.LHEHTCut)
    logger.info( "Adding upper LHE cut at %f", options.LHEHTCut )
    skimConds.append( "LHE_HTIncoming<%f"%options.LHEHTCut )

# final output directory 
storage_directory = os.path.join( options.targetDir, options.processingEra, str(options.year), options.skim, sample.name )
try:    #Avoid trouble with race conditions in multithreading
    os.makedirs(storage_directory)
    logger.info( "Created output directory %s.", storage_directory )
except:
    pass

## sort the list of files?
len_orig = len(sample.files)
sample = sample.split( n=options.nJobs, nSub=options.job)
logger.info( "fileBasedSplitting: Run over %i/%i files for job %i/%i."%(len(sample.files), len_orig, options.job, options.nJobs))
logger.debug("fileBasedSplitting: Files to be run over:\n%s", "\n".join(sample.files) )

# turn on all branches to be flexible for filter cut in skimCond etc.
sample.chain.SetBranchStatus("*",1)

# B tagging SF
b_tagger = "DeepCSV"
from Analysis.Tools.BTagEfficiency import BTagEfficiency
btagEff = BTagEfficiency( fastSim = False, year=options.year, tagger=b_tagger )

# tmp_output_directory
tmp_output_directory  = os.path.join( user.postprocessing_tmp_directory, "%s_%i_%s_%s_%s"%(options.processingEra, options.year, options.skim, sample.name, str(uuid.uuid3(uuid.NAMESPACE_OID, sample.name))))  

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

target_outfilename = os.path.join(storage_directory, sample.name + '.root') 
filename, ext      = os.path.splitext( os.path.join(tmp_output_directory, sample.name + '.root') )
outfilename        = filename+ext

if not options.overwrite:
    if os.path.isfile(target_outfilename):
        logger.info( "Output file %s found.", target_outfilename)
        if checkRootFile( target_outfilename, checkForObjects=["Events"] ) and deepCheckRootFile( target_outfilename ) and deepCheckWeight( target_outfilename ):
            logger.info( "File already processed. Source: File check ok! Skipping." ) # Everything is fine, no overwriting
            sys.exit(0)
        else:
            logger.info( "File corrupt. Removing file from target." )
            os.remove( target_outfilename )
            logger.info( "Reprocessing." )
    else:
        logger.info( "Sample not processed yet." )
        logger.info( "Processing." )
else:
    logger.info( "Overwriting.")

# relocate original
sample.copy_files( os.path.join(tmp_output_directory, "input") )

# Decision based on sample name -> whether TTJets or TTLep is in the sample name
isTT = sample.name.startswith("TTJets") or sample.name.startswith("TTLep") or sample.name.startswith("TT_pow")

if sample.name.startswith("TTLep"):
    sample.topScaleF = 1.002 ## found to be universal for years 2016-2018, and in principle negligible

if options.doCRReweighting:
    from tWZ.Tools.colorReconnectionReweighting import getCRWeight, getCRDrawString
    logger.info( "Sample will have CR reweighting." )
    selectionString = "&&".join(skimConds)
    #norm = sample.getYieldFromDraw( selectionString = selectionString, weightString = "genWeight" )
    norm = float(sample.chain.GetEntries(selectionString))
    CRScaleF = sample.getYieldFromDraw( selectionString = selectionString, weightString = getCRDrawString() )
    CRScaleF = CRScaleF['val']/norm#['val']
    logger.info(" Using norm: %s"%norm )
    logger.info(" Found CRScaleF: %s"%CRScaleF)
else:
    CRScaleF = 1
    logger.info( "Sample will NOT have CR reweighting. CRScaleF=%f",CRScaleF )

#branches to be kept for data and MC
branchKeepStrings_DATAMC = [\
    "run", "luminosityBlock", "event", "fixedGridRhoFastjetAll", "PV_npvs", "PV_npvsGood",
    "MET_*",
    "CaloMET_*",
    "RawMET_phi", "RawMET_pt", "RawMET_sumEt",
    "Flag_*",
    "nJet", "Jet_*",
    "nElectron", "Electron_*",
    "nMuon", "Muon_*",
]
branchKeepStrings_DATAMC += ["HLT_*", "PV_*"]

if options.year == 2017:
    branchKeepStrings_DATAMC += [ "METFixEE2017_*" ]

#branches to be kept for MC samples only
branchKeepStrings_MC = [ "Generator_*",  "genWeight", "Pileup_nTrueInt","GenMET_*", "nISR", "nGenJet", "GenJet_*"]
#branchKeepStrings_MC.extend([ "*LHEScaleWeight", "*LHEPdfWeight", "LHEWeight_originalXWGTUP"]) # don't exist in QCD MC?

#branches to be kept for data only
branchKeepStrings_DATA = [ ]

# Jet variables to be read from chain
jetCorrInfo = []
jetMCInfo   = ['genJetIdx/I','hadronFlavour/I']

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

jetVars        = ['pt/F', 'chEmEF/F', 'chHEF/F', 'neEmEF/F', 'neHEF/F', 'rawFactor/F', 'eta/F', 'phi/F', 'jetId/I', 'btagDeepB/F', 'btagDeepFlavB/F', 'btagCSVV2/F', 'area/F', 'pt_nom/F', 'corr_JER/F'] + jetCorrInfo
if isMC:
    jetVars    += jetMCInfo
jetVarNames    = [x.split('/')[0] for x in jetVars]
genLepVars     = ['pt/F', 'phi/F', 'eta/F', 'pdgId/I', 'genPartIdxMother/I', 'status/I', 'statusFlags/I'] # some might have different types
genLepVarNames = [x.split('/')[0] for x in genLepVars]

muVars         = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F', 'pfRelIso03_all/F', 'pfRelIso04_all/F', 'mvaTOP/F', 'sip3d/F','convVeto/I','dxy/F','dz/F','charge/I','mediumId/I','index/I', 'mT/F', 'hybridIso/F']
eleVars        = ['pt/F','eta/F','phi/F','pdgId/I','cutBased/I','miniPFRelIso_all/F', 'pfRelIso03_all/F', 'mvaTOP/F', 'sip3d/F','convVeto/I','dxy/F','dz/F','charge/I','deltaEtaSC/F','mediumId/I','index/I', 'mT/F', 'hybridIso/F', 'sieie/F','hoe/F','eInvMinusPInv/F', 'mvaFall17V2noIso_WP80/O']

muVarNames     = [x.split('/')[0] for x in muVars] 
eleVarNames    = [x.split('/')[0] for x in eleVars]

read_variables = map(TreeVariable.fromString, [ 'MET_pt/F', 'MET_phi/F', 'run/I', 'luminosityBlock/I', 'event/l', 'PV_npvs/I', 'PV_npvsGood/I'] )
if options.year == 2017:
    read_variables += map(TreeVariable.fromString, [ 'METFixEE2017_pt/F', 'METFixEE2017_phi/F', 'METFixEE2017_pt_nom/F', 'METFixEE2017_phi_nom/F'])
else:
    read_variables += map(TreeVariable.fromString, [ 'MET_pt_nom/F', 'MET_phi_nom/F' ])

read_variables += [ TreeVariable.fromString('nPhoton/I'),
                    VectorTreeVariable.fromString('Photon[pt/F,eta/F,phi/F,mass/F,cutBased/I,pdgId/I]') if (options.year == 2016) else VectorTreeVariable.fromString('Photon[pt/F,eta/F,phi/F,mass/F,cutBasedBitmap/I,pdgId/I]') ]

new_variables = [ 'weight/F', 'triggerDecision/I', 'year/I']
if isMC:
    read_variables += [ TreeVariable.fromString('Pileup_nTrueInt/F') ]
    # reading gen particles for top pt reweighting
    read_variables.append( TreeVariable.fromString('nGenPart/I') )
    read_variables.append( VectorTreeVariable.fromString('GenPart[pt/F,mass/F,phi/F,eta/F,pdgId/I,genPartIdxMother/I,status/I,statusFlags/I]', nMax=200 )) # default nMax is 100, which would lead to corrupt values in this case
    read_variables.append( TreeVariable.fromString('genWeight/F') )
    read_variables.append( TreeVariable.fromString('nGenJet/I') )
    read_variables.append( VectorTreeVariable.fromString('GenJet[pt/F,eta/F,phi/F]' ) )
    new_variables.extend([ 'reweight_nISR/F', 'reweight_nISRUp/F', 'reweight_nISRDown/F', 'reweightPU/F','reweightPUUp/F','reweightPUDown/F', 'reweightPUVUp/F','reweightPUVVUp/F', 'reweightPUVDown/F', 'reweightL1Prefire/F', 'reweightL1PrefireUp/F', 'reweightL1PrefireDown/F'])
#    if not options.skipGenLepMatching:
#        new_variables.append( TreeVariable.fromString( 'nGenLep/I' ) )
#        new_variables.append( 'GenLep[%s]'% ( ','.join(genLepVars) ) )
    if options.doCRReweighting:
        new_variables.append('reweightCR/F')

read_variables += [\
    TreeVariable.fromString('nElectron/I'),
    VectorTreeVariable.fromString('Electron[pt/F,eta/F,phi/F,pdgId/I,cutBased/I,miniPFRelIso_all/F,pfRelIso03_all/F,sip3d/F,lostHits/b,mvaFall17V2Iso_WP80/O,mvaFall17V2Iso_WP90/O,convVeto/O,dxy/F,dz/F,charge/I,deltaEtaSC/F,vidNestedWPBitmap/I,mvaTOP/F,jetIdx/I,jetRelIso/F,sieie/F,hoe/F,eInvMinusPInv/F,mvaFall17V2noIso_WP80/O]'),
    TreeVariable.fromString('nMuon/I'),
    VectorTreeVariable.fromString('Muon[pt/F,eta/F,phi/F,pdgId/I,mediumId/O,miniPFRelIso_all/F,pfRelIso03_all/F,pfRelIso04_all/F,sip3d/F,dxy/F,dz/F,charge/I,mvaTOP/F,looseId/O,jetIdx/I,jetRelIso/F,mvaId/b]'),
    TreeVariable.fromString('nJet/I'),
    VectorTreeVariable.fromString('Jet[%s]'% ( ','.join(jetVars) ) ) ]


read_variables.extend(
    [
        TreeVariable.fromString("{trigger}/O".format(trigger =t)) for t in options.trigger
    ]
)


if sample.isData: new_variables.extend( ['jsonPassed/I','isData/I'] )
new_variables.extend( ['nBTag/I'] )

# ids
mu_ids  = ['FOmvaTOPT', 'mvaTOPT', 'hybridIso', 'looseHybridIso']
ele_ids = ['FOmvaTOPT', 'mvaTOPT', 'hybridIso', 'looseHybridIso']

L_T_pairs = {'looseHybridIso':'hybridIso', 'FOmvaTOPT':'mvaTOPT'}

for mu_id in mu_ids:
    new_variables.append( 'mu_%s[%s]'% (mu_id, ','.join(muVars) + ',mvaTOPWP/I,pt_corr/F,deepJet/F' ) )
    new_variables.append( 'nmu_%s/I'%mu_id )
for ele_id in ele_ids:
    new_variables.append( 'ele_%s[%s]'% ( ele_id, ','.join(eleVars) + ',mvaTOPWP/I,pt_corr/F,deepJet/F,lostHitsI/I' ) )
    new_variables.append( 'nele_%s/I'%ele_id )

cleaning_vars  = ["cleaned_mu_%s/I"%mu_id for mu_id in mu_ids]
cleaning_vars += ["cleaned_ele_%s/I"%ele_id for ele_id in ele_ids]
new_variables += [\
    'overlapRemoval/I',
    'JetGood[%s]'% ( ','.join(jetVars+['index/I']+cleaning_vars) + ',genPt/F' ),
    'met_pt/F', 'met_phi/F', 'met_pt_min/F']

if options.checkTTGJetsOverlap:
    new_variables.extend( ['TTGJetsEventType/I'] )

# Btag weights Method 1a
if isMC:
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
    #from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetUncertainties   import jetmetUncertaintiesProducer
    #from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetRecalib            import jetRecalib
    from PhysicsTools.NanoAODTools.postprocessing.modules.common.ISRcounter        import ISRcounter
    
    logger.info("Preparing nanoAOD postprocessing")
    logger.info("Will put files into directory %s", tmp_output_directory)
    cut = '&&'.join(skimConds)
    # year specific JECs 
    if options.year == 2016:
        JER                 = "Summer16_25nsV1_MC"          if not sample.isData else "Summer16_25nsV1_DATA"
        JERera              = "Summer16_25nsV1"

    elif options.year == 2017:
        JER                 = "Fall17_V3_MC"                if not sample.isData else "Fall17_V3_DATA"
        JERera              = "Fall17_V3"

    elif options.year == 2018:
        JER                 = "Autumn18_V1_MC"              if not sample.isData else "Autumn18_V1_DATA"
        JERera              = "Autumn18_V1"

    if options.overwriteJEC is not None:
        JEC = options.overwriteJEC

    logger.info("Using JERs for MET significance: %s", JER)
    
    from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jetmetHelperRun2 import *
    METBranchName = 'MET' if not options.year == 2017 else 'METFixEE2017'

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
        
        if not sample.isData:
            modules.append( ISRcounter() )

        # need a hash to avoid data loss
        file_hash = str(hash(f))
        p = PostProcessor(tmp_output_directory, [f], cut=cut, modules=modules, postfix="_for_%s_%s"%(sample.name, file_hash))
        if not options.reuseNanoAOD:
            p.run()
        newFileList += [tmp_output_directory + '/' + f.split('/')[-1].replace('.root', '_for_%s_%s.root'%(sample.name, file_hash))]
    logger.info("Done. Replacing input files for further processing.")
    
    sample.files = newFileList
    sample.clear()

# Define a reader
reader = sample.treeReader( \
    variables = read_variables,
    selectionString = "&&".join(skimConds)
    )

# using miniRelIso 0.2 as baseline 

eleSelector_ = { ele_id:eleSelector( ele_id, year = options.year ) for ele_id in ele_ids }
muSelector_  = { mu_id:muonSelector( mu_id,  year = options.year ) for mu_id in mu_ids }

mothers = {"D":0, "B":0}
grannies_D = {}
grannies_B = {}

def filler( event ):
    # shortcut
    r = reader.event
    #workaround  = (r.run, r.luminosityBlock, r.event) # some fastsim files seem to have issues, apparently solved by this.
    event.isData = s.isData
    event.year   = options.year
    event.overlapRemoval = 1 
        
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

    if isMC and hasattr(r, "Pileup_nTrueInt"):
        event.reweightPU     = nTrueInt_puRW       ( r.Pileup_nTrueInt ) # is this correct?
        event.reweightPUDown = nTrueInt_puRWDown   ( r.Pileup_nTrueInt )
        event.reweightPUVDown= nTrueInt_puRWVDown  ( r.Pileup_nTrueInt )
        event.reweightPUUp   = nTrueInt_puRWUp     ( r.Pileup_nTrueInt )
        event.reweightPUVUp  = nTrueInt_puRWVUp    ( r.Pileup_nTrueInt )
        event.reweightPUVVUp = nTrueInt_puRWVVUp   ( r.Pileup_nTrueInt )

    # Trigger Decision
    event.triggerDecision = int(treeFormulas['triggerDecision']['TTreeFormula'].EvalInstance())

    allSlimmedJets      = getJets(r)
    allSlimmedPhotons   = getPhotons(r, year=options.year)
    if options.year == 2018:
        event.reweightL1Prefire, event.reweightL1PrefireUp, event.reweightL1PrefireDown = 1., 1., 1.
    else:
        event.reweightL1Prefire, event.reweightL1PrefireUp, event.reweightL1PrefireDown = L1PW.getWeight(allSlimmedPhotons, allSlimmedJets)

    # store the correct MET (EE Fix for 2017, MET_min as backup in 2017)
    
    if options.year == 2017:# and not options.fastSim:
        # v2 recipe. Could also use our own recipe
        event.met_pt    = r.METFixEE2017_pt_nom
        event.met_phi   = r.METFixEE2017_phi_nom
        #event.met_pt_min = r.MET_pt_min not done anymore
    else:
        event.met_pt    = r.MET_pt_nom 
        event.met_phi   = r.MET_phi_nom

        event.met_pt_min = 0

    all_jets     = getJets(r, jetVars=jetVarNames)

    all_muons      = getGoodMuons    (r, mu_selector = None )
    for m in all_muons:
        if m['jetIdx']>=0:
            m['deepJet'] = all_jets[m['jetIdx']]['btagDeepFlavB']
        else:
            m['deepJet'] = -1#float('nan') 

    for mu_id in mu_ids:
        muons      = filter( muSelector_[mu_id], all_muons)
        for m in muons:
            m['pdgId']      = int( -13*m['charge'] )
            m['mT']         = sqrt( 2*event.met_pt*m['pt']*(1-cos(event.met_phi-m['phi'])) )
            m['hybridIso']  = m['pfRelIso04_all']*min(25, m['pt'])
            m['mvaTOPWP']   = mvaTopWP(m['mvaTOP'], m['pdgId'])
            if L_T_pairs.has_key(mu_id): # this mu_id has a tight ID that corresponds to it
                is_tight = muSelector_[L_T_pairs[mu_id]](m)
                if is_tight:# or m['jetIdx']<0:
                    m['pt_corr'] = m['pt'] 
                else:
                    m['pt_corr'] = 0.8*m['pt']*(1+m['jetRelIso'])
            else:
                    m['pt_corr'] = m['pt'] 
        
        fill_vector_collection( event, "mu_%s"%mu_id,  muVarNames,  muons)
        setattr(event, "nmu_%s"%mu_id, len(muons) )

        clean_jets,_ = cleanJetsAndLeptons( all_jets, muons ) 
        for jet in all_jets:
            jet["cleaned_mu_%s"%mu_id] = (jet in clean_jets)

    all_electrons = getGoodElectrons(r, ele_selector = None)
    #print all_electrons
    for e in all_electrons:
        e['lostHitsI'] = ord(e['lostHits'])
        if e['jetIdx']>=0:
            e['deepJet'] = all_jets[e['jetIdx']]['btagDeepFlavB']
        else:
            e['deepJet'] = -1#float('nan') 
    for ele_id in ele_ids:
        #print ele_id, len(all_electrons)
        electrons = filter(eleSelector_[ele_id], all_electrons)
        #print ele_id, len(electrons)
        for e in electrons:
            e['pdgId']      = int( -11*e['charge'] )
            #e['mT']         = sqrt( 2*event.met_pt*e['pt']*(1-cos(event.met_phi-e['phi'])) )
            e['hybridIso']  = e['pfRelIso03_all']*min(25, e['pt'])
            e['mvaTOPWP']   = mvaTopWP(e['mvaTOP'], e['pdgId'])
            if L_T_pairs.has_key(ele_id): # this ele_id has a tight ID that corresponds to it
                is_tight = eleSelector_[L_T_pairs[ele_id]](e)
                if is_tight:# or e['jetIdx']<0:
                    e['pt_corr'] = e['pt'] 
                else:
                    e['pt_corr'] = 0.8*e['pt']*(1+e['jetRelIso'])
            else:
                    e['pt_corr'] = e['pt'] 
                

        fill_vector_collection( event, "ele_%s"%ele_id, eleVarNames, electrons)
        setattr(event, "nele_%s"%ele_id, len(electrons) )

        clean_jets,_ = cleanJetsAndLeptons( all_jets, electrons ) 
        for jet in all_jets:
            jet["cleaned_ele_%s"%ele_id] = (jet in clean_jets)

    clean_jets_acc = filter(lambda j:abs(j['eta'])<2.4, all_jets)
    jets         = filter(lambda j:j['pt']>30, clean_jets_acc)
    bJets        = filter(lambda j:isBJet(j, tagger=b_tagger, year=options.year) and abs(j['eta'])<=2.4    , jets)
    nonBJets     = filter(lambda j:not ( isBJet(j, tagger=b_tagger, year=options.year) and abs(j['eta'])<=2.4 ), jets)

    # Filling jets
    maxNJet = 100
    store_jets = jets #if not options.keepAllJets else soft_jets + jets
    store_jets = store_jets[:maxNJet]
    store_jets.sort( key = lambda j:-j['pt'])
    event.nJetGood   = len(store_jets)
    for iJet, jet in enumerate(store_jets):
        event.JetGood_index[iJet] = jet['index']
        if isMC:
            if store_jets[iJet]['genJetIdx'] >= 0:
                if r.nGenJet<maxNJet:
                    try:
                        event.JetGood_genPt[iJet] = r.GenJet_pt[store_jets[iJet]['genJetIdx']]
                    except IndexError:
                        event.JetGood_genPt[iJet] = -1
                else:
                    event.JetGood_genPt[iJet] = -1
        for b in jetVarNames + [x.split('/')[0] for x in cleaning_vars]:
            getattr(event, "JetGood_"+b)[iJet] = jet[b]
        #getattr(event, "JetGood_pt")[iJet] = jet['pt']

    if isMC and options.doCRReweighting:
        event.reweightCR = getCRWeight(event.nJetGood)

    event.nBTag      = len(bJets)

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
    nEvents = clonedTree.GetEntries()
    if nEvents==0:
        logger.info( "No Events found. Skip." )
        continue
        
    clonedEvents += nEvents

    # Add the TTreeFormulas
    for formula in treeFormulas.keys():
        treeFormulas[formula]['TTreeFormula'] = ROOT.TTreeFormula(formula, treeFormulas[formula]['string'], clonedTree )

    # Clone the empty maker in order to avoid recompilation at every loop iteration
    maker = treeMaker_parent.cloneWithoutCompile( externalTree = clonedTree )

    maker.start()
    # Do the thing
    reader.start()

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

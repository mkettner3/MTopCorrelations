# Standard
import os
import ROOT
import subprocess

# RootTools
from RootTools.core.standard                     import *

# Analysis
from Analysis.Tools.helpers import deepCheckRootFile, checkRootFile

# User specific
from tWZ.Tools.user import postprocessing_output_directory

def get_parser():
    ''' Argument parser for post-processing module.
    '''
    import argparse
    argParser = argparse.ArgumentParser(description = "Argument parser for nanoPostProcessing")
    argParser.add_argument('--file',       action='store', type=str, default='nanoPostProcessing_Summer16', help="postprocessing sh file to check")
    argParser.add_argument('--createExec', action='store_true', help="create .sh file with missing files?")
    argParser.add_argument('--overwrite',  action='store_true', help="overwrite existing missingFiles.sh file?")
    argParser.add_argument('--check',      action='store', nargs='?',  default='normal', choices=['normal', 'deep', 'None'], help="Perform check?")
    argParser.add_argument('--logLevel',   action='store', nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
    argParser.add_argument('--data_directory',   action='store', nargs='?',  default = postprocessing_output_directory, help="Which directory to check?" )
    return argParser

args = get_parser().parse_args()

# Logging
if __name__=="__main__":
    import Analysis.Tools.logger as logger
    logger = logger.get_logger(args.logLevel, logFile = None )
    import RootTools.core.logger as logger_rt
    logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

else:
    import logging
    logger = logging.getLogger(__name__)

if args.file.endswith(".sh"):
    args.file = args.file.rstrip(".sh")

def filterEmpty( strList ):
    return list( filter ( bool, strList ) )

def getDataDictList( filepath ):
    ''' Read postprocessing sh file and format it to dictionary
    '''
    with open( filepath, 'r' ) as f:
        ppLines = f.readlines()

    ppLines = [ line for line in ppLines if line.startswith('python') ]

    dictList = []
    for line in ppLines:
        skim    = filterEmpty( line.split("--skim ")[1].split(" ") )[0]
        year    = filterEmpty( line.split("--year ")[1].split(" ") )[0]
        dir     = filterEmpty( line.split("--processingEra ")[1].split(" ") )[0]
        sample  = filterEmpty( line.split("--sample ")[1].split(" ") )[0]
        command = line
        # find sample arguments
        try:
            sample_args = line.split("--sample")[1].split('--')[0].split('#')[0].rstrip().lstrip().split()
        except IndexError as e:
            logger.error("No sample argument in line? Reading: %s" % line)
            raise e
        if len(sample_args)>1:
            sample += "_comb"
        if "#SPLIT" not in line: 
            nFiles=1
        else:
            nFiles = filterEmpty( line.split("#SPLIT")[1].split(" ") )[0].split("\n")[0]
        dictList.append( { "skim":skim, "year":int(year), "dir":dir, "sample":sample, "nFiles":int(nFiles), "command":command} )

    return dictList

# Load File
logger.info( "Now running on pp file %s and checking in directory %s", args.file, args.data_directory)
file          = os.path.expandvars( "$CMSSW_BASE/src/tWZ/postProcessing/%s.sh" % args.file )
dictList      = getDataDictList( file )
isData        = "Run" in args.file
execCommand   = []
for ppEntry in dictList:
    sample = ppEntry['sample']
    logger.debug("Checking sample %s" %sample)

    postfix = ""
    # check whether we do reduction:
    if "--reduceSizeBy" in ppEntry["command"]:
        postfix += "_redBy%i" % int(ppEntry["command"].split('--reduceSizeBy')[1].lstrip().split(' ')[0]) 
    dirPath = os.path.join( args.data_directory, ppEntry["dir"], ppEntry["skim"], sample+postfix  )
    # check whether we have an LHE HT cut:
    if "--LHEHTCut" in ppEntry["command"]:
        postfix += "_lheHT%i" % int(ppEntry["command"].split('--LHEHTCut')[1].lstrip().split(' ')[0]) 
    dirPath = os.path.join( args.data_directory, ppEntry["dir"], str(ppEntry["year"]), ppEntry["skim"], sample+postfix  )

    # find all root files in subdirectory

    allFiles = os.listdir(dirPath) if os.path.exists( dirPath ) else []
    allRootFiles = filter( lambda f: f.endswith('.root') and not f.startswith('nanoAOD'), allFiles )
    nanoAODFiles = filter( lambda f: f.endswith('.root') and f.startswith('nanoAOD'), allFiles )
    prefix = ''

    nanoAOD_list = []
    if len( nanoAODFiles )>0:
        logger.warning( "Found nanoAODFiles. Shouldn't be here: %s", ",".join( nanoAODFiles ) )
        try:
            nanoAOD_list = list(set([ int(f.rstrip('.root').split('_')[-2]) for f in nanoAODFiles ]))
        except ValueError: #possibly 'SPLIT1'
            nanoAOD_list = [0]
    
    if len(allRootFiles)>0:
        rootFiles  = [ os.path.join(dirPath, filename) for filename in allRootFiles]

        if args.check!='None':
            if args.check=='deep':
                check_f = lambda f: checkRootFile(prefix+f, ["Events"]) and deepCheckRootFile( prefix+f ) 
            elif args.check=='normal':
                check_f = lambda f: checkRootFile(prefix+f, ["Events"])

            len_before = len(rootFiles)
            rootFiles = filter(check_f, rootFiles )
            if len_before>len(rootFiles):
                logger.warning( "Path: %s : %i/%i files failed %s check", dirPath, len_before-len(rootFiles), len_before, args.check)
            elif len_before==len(rootFiles):
                logger.debug( "Sample %s: All %i files passed %s check", sample, len(rootFiles), args.check)
    else: 
        rootFiles = []
        logger.warning( "Path does not exist or no files found: %s", dirPath )


    if not rootFiles:
        logger.info("Sample %s not processed" %sample)
        if args.createExec:
            execCommand += [ ppEntry["command"] ]
        continue

    # more than one files
    if ppEntry['nFiles']>1:
        try:
            job_ids = [ int(item.rstrip('.root').split('_')[-1]) for item in rootFiles ]
        except:
            print rootFiles
    # one file
    elif ppEntry['nFiles']==1 and rootFiles[0].endswith(sample+'.root'):
        job_ids = [0]
    
    missing_ids = [ id for id in range(ppEntry['nFiles']) if ( id not in job_ids or id in nanoAOD_list ) ]
    if args.createExec:
        for missing_id in missing_ids:
            execCommand += [ ppEntry["command"].split("#SPLIT")[0] + "--nJobs %s --job %s"%(ppEntry["nFiles"], missing_id) ]

    if len(missing_ids)>0:
        logger.info( "Sample %s: Found %i missing IDs: %s", sample, len(missing_ids), ",".join( map(str, missing_ids)) )
    else:
        logger.info( "Sample %s seems OK. %i jobs", sample, ppEntry['nFiles'])
if args.createExec:
    with open( "missingFiles.sh", 'w' if args.overwrite else "a" ) as f:
        for line in execCommand:
            f.write(line.split("\n")[0] + "\n")
        f.write("\n")


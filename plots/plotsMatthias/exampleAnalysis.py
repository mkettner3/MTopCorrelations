#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('delete.png')
import itertools
import copy
import array
import operator
from math                                import sqrt, cos, sin, pi, atan2, cosh, exp

# RootTools
from RootTools.core.standard             import *

# MTopCorrelations
from MTopCorrelations.Tools.user                      import plot_directory
from MTopCorrelations.Tools.cutInterpreter            import cutInterpreter
from MTopCorrelations.Tools.objectSelection           import cbEleIdFlagGetter, vidNestedWPBitMapNamingList
from MTopCorrelations.Tools.objectSelection           import lepString
from MTopCorrelations.Tools.helpers          import getCollection

# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

import Analysis.Tools.syncer
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory', action='store', default='MTopCorrelations_v1')
argParser.add_argument('--selection',      action='store', default='nAK82p')
argParser.add_argument('--era',            action='store', type=str, default="UL2018")
args = argParser.parse_args()

################################################################################
# Logger
import MTopCorrelations.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

################################################################################
# Define the MC samples
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import *

mc = [UL2018.TTbar]
lumi_scale = 60

################################################################################
# Correlator Hist 
# hist = ROOT.TH1F("Correlator", "dR", 10, 0, 3)

################################################################################
# Text on the plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

################################################################################
# Functions needed specifically for this analysis routine

def drawObjects( plotData, lumi_scale ):
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if plotData else 'CMS Simulation'),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) '% ( lumi_scale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines]

def drawPlots(plots):
    for log in [False, True]:
        plot_directory_ = os.path.join(plot_directory, 'analysisPlots', args.plot_directory, args.era, ("_log" if log else ""), args.selection)
        for plot in plots:
            if not max(l.GetMaximum() for l in sum(plot.histos,[])): continue # Empty plot

            _drawObjects = []
            n_stacks=len(plot.histos)
            plotData=False
            if isinstance( plot, Plot):
                plotting.draw(plot,
                  plot_directory = plot_directory_,
                  ratio =  None,
                  logX = False, logY = log, sorting = True,
                  yRange = (0.03, "auto") if log else (0.001, "auto"),
                  scaling = {},
                  legend = ( (0.18,0.88-0.03*sum(map(len, plot.histos)),0.9,0.88), 2),
                  drawObjects = drawObjects( plotData , lumi_scale ) + _drawObjects,
                  copyIndexPHP = True, extensions = ["png", "pdf", "root"],
                )

################################################################################
# Define sequences
sequence       = []


################################################################################
# Read variables

read_variables = [
    "nGenPart/I",
    "GenPart[pt/F,eta/F,phi/F,m/F,pdgId/I,mompdgId/i]",
    # "nGenJetAK8/I",
    # "GenJetAK8[pt/F,eta/F,phi/F,mass/F]",
    "nPFJetAK8/I",
    "PFJetAK8[pt/F,eta/F,phi/F,mass/F]",
    "nPFJetAK8_cons/I",
    VectorTreeVariable.fromString( "PFJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000),
    
]

################################################################################
# Set up plotting
# weight_ = lambda event, sample: event.weight if sample.isData else event.weight*lumi_scale/1000.
# weight_ = lambda event, sample: 1. if sample.isData else lumi_scale/1000.
weight_ = lambda event, sample: 1.

for sample in mc: 
    sample.style = styles.fillStyle(sample.color)

for sample in mc:
    sample.weight = lambda event, sample: 1.
    
stack = Stack(mc)

# Use some defaults
Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

################################################################################
# Now define the plots

plots = []

plots.append(Plot(
    name = "nAK8",
    texX = 'Number of AK8 jets', texY = 'Number of Events',
    attribute = lambda event, sample: event.nPFJetAK8,
    binning=[11, -0.5, 10.5],
))

plots.append(Plot(
    name = "mjet",
    texX = 'm_{jet} [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.PFJetAK8_mass[0],
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "ptjet",
    texX = 'Leading AK8 jet p_{T} [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.PFJetAK8_pt[0],
    binning=[25, 0., 500.],
))


plotting.fill(plots, read_variables = read_variables, sequence = sequence)

drawPlots(plots)
    
# print hist.Integral()    
    
logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

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
# import itertools
# import copy
# import array
# import operator
# from math                                import sqrt, cos, sin, pi, atan2, cosh, exp

# RootTools
from RootTools.core.standard             import *

# MTopCorrelations
from MTopCorrelations.Tools.user                      import plot_directory
from MTopCorrelations.Tools.cutInterpreter            import cutInterpreter
from MTopCorrelations.Tools.objectSelection           import cbEleIdFlagGetter, vidNestedWPBitMapNamingList
from MTopCorrelations.Tools.objectSelection           import lepString
from MTopCorrelations.Tools.helpers          import getCollection

# Analysis
# from Analysis.Tools.helpers              import deltaPhi, deltaR
# from Analysis.Tools.puProfileCache       import *
# from Analysis.Tools.puReweighting        import getReweightingFunction
# from Analysis.Tools.leptonJetArbitration     import cleanJetsAndLeptons

# import Analysis.Tools.syncer
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory', action='store', default='MTopCorrelations_v1')
argParser.add_argument('--selection',      action='store', default='nAK81p-AK8pt')
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
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018

mc = [UL2018.TTbar]
lumi_scale = 60

################################################################################
# Correlator Hist
hist = ROOT.TH1F("Correlator", "dR", 10, 0, 3)

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


def getJetConstituents(event, idx):
    constituents = []
    for i in range(event.nGenJetAK8_cons):
        if event.GenJetAK8_cons_jetIndex[i] == idx:
            part = ROOT.TLorentzVector()
            part.SetPtEtaPhiM(event.GenJetAK8_cons_pt[i],event.GenJetAK8_cons_eta[i],event.GenJetAK8_cons_phi[i],event.GenJetAK8_cons_mass[i])
            constituents.append(part)
    def getPt(part):
        return part.Pt()
    constituents.sort(key=getPt, reverse=True)
    return constituents


def gen_tops(event, sample):
    top_vec = ROOT.TLorentzVector()
    anti_top_vec = ROOT.TLorentzVector()
    for i in range(event.nGenPart):
        if event.GenPart_pdgId[i] == 6:
            top_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
        elif event.GenPart_pdgId[i] == -6:
            anti_top_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])

    top_lep = ROOT.TLorentzVector()
    top_had = ROOT.TLorentzVector()
    for i in range(event.nGenPart):
         if abs(event.GenPart_pdgId[i]) in [11, 13, 15]:
            if event.GenPart_grmompdgId[i] == 6:
                top_lep = top_vec
                top_had = anti_top_vec
            elif event.GenPart_grmompdgId[i] == -6:
                top_lep = anti_top_vec
                top_had = top_vec

    # event.top_had_pt = top_had.Pt()

    quark_vec = ROOT.TLorentzVector()
    anti_quark_vec = ROOT.TLorentzVector()
    bottom_vec = ROOT.TLorentzVector()
    for i in range(event.nGenPart):
        if event.GenPart_pdgId[i] in range(1, 7) and abs(event.GenPart_grmompdgId[i]) == 6:
            quark_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
        elif event.GenPart_pdgId[i] in range(-1, -7, -1) and abs(event.GenPart_grmompdgId[i]) == 6:
            anti_quark_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
        elif event.GenPart_pdgId[i] == 5 and abs(event.GenPart_mompdgId[i]) == 6:
            bottom_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])

    jets = []
    for i in range(event.nGenJetAK8):
        jets.append(ROOT.TLorentzVector())
        jets[-1].SetPtEtaPhiM(event.GenJetAK8_pt[i], event.GenJetAK8_eta[i], event.GenJetAK8_phi[i], event.GenJetAK8_mass[i])
    deltas_t = [jets[j].DeltaR(top_had) for j in range(event.nGenJetAK8)]
    event.delta_min = min(deltas_t)
    nearest_jet_idx_t = np.argmin(deltas_t)
    event.nearest_jet_mass_t = jets[nearest_jet_idx_t].M()
    event.nearest_jet_pt = jets[nearest_jet_idx_t].Pt()
    event.nearest_jet_eta = jets[nearest_jet_idx_t].Eta()
    event.nearest_jet_phi = jets[nearest_jet_idx_t].Phi()

    delta_q = jets[nearest_jet_idx_t].DeltaR(quark_vec)
    delta_aq = jets[nearest_jet_idx_t].DeltaR(anti_quark_vec)
    delta_b = jets[nearest_jet_idx_t].DeltaR(bottom_vec)

    event.merged_top = False
    event.merged_q_aq_b = False

    if event.delta_min < 0.8:
        event.merged_top = True

        if delta_q < 0.8 and delta_aq < 0.8 and delta_b < 0.8:
            event.merged_q_aq_b = True

    event.matched_jet_nCons = 0
    event.matched_jet_Cons_pt = []
    for k in range(event.nGenJetAK8_cons):
        if event.GenJetAK8_cons_jetIndex[k] == nearest_jet_idx_t:
            event.matched_jet_nCons += 1
            event.matched_jet_Cons_pt.append(event.GenJetAK8_cons_pt[k])

    event.matched_jet_Cons_pt_ratio = [elem/event.nearest_jet_pt for elem in event.matched_jet_Cons_pt]

    matched_jet_cons = getJetConstituents(event=event, idx=nearest_jet_idx_t)

    delta_delta = 0.02
    triplet = [ROOT.TLorentzVector()]*3
    numb_of_particles = 1000
    if len(matched_jet_cons) >= numb_of_particles:
        for i in range(numb_of_particles):
            for j in range(i+1, numb_of_particles):
                for k in range(j+1, numb_of_particles):

                    triplet[0] = matched_jet_cons[i]
                    triplet[1] = matched_jet_cons[j]
                    triplet[2] = matched_jet_cons[k]
                    delta_0_1 = triplet[0].DeltaR(triplet[1])
                    delta_0_2 = triplet[0].DeltaR(triplet[2])
                    delta_1_2 = triplet[1].DeltaR(triplet[2])
                    if abs(delta_0_1-delta_0_2) < delta_delta:
                        if abs(delta_1_2-delta_0_1) < delta_delta:
                            if abs(delta_0_2-delta_1_2) < delta_delta:
                                w = (event.GenJetAK8_cons_pt[i]*event.GenJetAK8_cons_pt[j]*event.GenJetAK8_cons_pt[k])
                                w = w**2
                                w = w / (event.nearest_jet_pt**(3*2))

                                zeta = (delta_0_1 + delta_0_2 + delta_1_2) / 3

                                hist.Fill(zeta, w)


sequence.append(gen_tops)

################################################################################
# Read variables

read_variables = [
    "nGenPart/I",
    "GenPart[pt/F,eta/F,phi/F,m/F,pdgId/I,mompdgId/I,grmompdgId/I]",
    "nGenJetAK8/I",
    "GenJetAK8[pt/F,eta/F,phi/F,mass/F]",
    "nPFJetAK8/I",
    "PFJetAK8[pt/F,eta/F,phi/F,mass/F]",
    "nPFJetAK8_cons/I",
    VectorTreeVariable.fromString( "PFJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000),
    "nGenJetAK8_cons/I",
    VectorTreeVariable.fromString( "GenJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000),
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
    name = "ptjet",
    texX = 'Leading AK8 jet p_{T} [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.PFJetAK8_pt[0],
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "delta_min",
    texX = 'Delta between nearest jet and hadronic top [rad]', texY = 'Number of Events',
    attribute = lambda event, sample: event.delta_min,
    binning=[25, 0., 7.],
))

plots.append(Plot(
    name = "jet_mass",
    texX = 'Mass of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_mass_t,
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "jet_pt",
    texX = 'P_{t} of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_pt,
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "jet_eta",
    texX = '#eta of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_eta,
    binning=[25, 3.5, 3.5],
))

plots.append(Plot(
    name = "jet_phi",
    texX = '#phi of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_phi,
    binning=[25, 3.5, 3.5],
))

plots.append(Plot(
    name = "cut_delta_min",
    texX = 'Delta between nearest jet and hadronic top [rad]', texY = 'Number of Events',
    attribute = lambda event, sample: event.delta_min if event.merged_top else float('nan'),
    binning=[25, 0., 7.],
))

plots.append(Plot(
    name = "cut_jet_mass",
    texX = 'Mass of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_mass_t if event.merged_top else float('nan'),
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "cut_jet_pt",
    texX = 'P_{t} of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_pt if event.merged_top else float('nan'),
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "cut_jet_eta",
    texX = '#eta of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_eta if event.merged_top else float('nan'),
    binning=[25, 3.5, 3.5],
))

plots.append(Plot(
    name = "cut_jet_phi",
    texX = '#phi of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_phi if event.merged_top else float('nan'),
    binning=[25, 3.5, 3.5],
))

plots.append(Plot(
    name = "cut_delta_min_2",
    texX = 'Delta between nearest jet and hadronic top [rad]', texY = 'Number of Events',
    attribute = lambda event, sample: event.delta_min if event.merged_q_aq_b else float('nan'),
    binning=[25, 0., 7.],
))

plots.append(Plot(
    name = "cut_jet_mass_2",
    texX = 'Mass of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_mass_t if event.merged_q_aq_b else float('nan'),
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "cut_jet_pt_2",
    texX = 'P_{t} of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_pt if event.merged_q_aq_b else float('nan'),
    binning=[25, 0., 500.],
))

plots.append(Plot(
    name = "cut_jet_eta_2",
    texX = '#eta of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_eta if event.merged_q_aq_b else float('nan'),
    binning=[25, 3.5, 3.5],
))

plots.append(Plot(
    name = "cut_jet_phi_2",
    texX = '#phi of jet next to hadronic top [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.nearest_jet_phi if event.merged_q_aq_b else float('nan'),
    binning=[25, 3.5, 3.5],
))

plots.append(Plot(
    name = "number_of_cons",
    texX = 'Number of constituents in matched jet', texY = 'Number of Events',
    attribute = lambda event, sample: event.matched_jet_nCons,
    binning=[25, 0., 150.],
))

plots.append(Plot(
    name = "cons_pt",
    texX = 'P_{t} of constituents of matched jet [GeV]', texY = 'Number of Events',
    attribute = lambda event, sample: event.matched_jet_Cons_pt,
    binning=[25, 0., 50.],
))

plots.append(Plot(
    name = "cons_pt_ratio",
    texX = 'P_{t} ratio of constituents of matched jet', texY = 'Number of Events',
    attribute = lambda event, sample: event.matched_jet_Cons_pt_ratio,
    binning=[25, 0., 0.25],
))


plotting.fill(plots, read_variables = read_variables, sequence = sequence)

drawPlots(plots)

# print hist.Integral()

f = ROOT.TFile('correlator.root', 'RECREATE')
f.cd()
hist.Write('correlator_hist')
f.Close()

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

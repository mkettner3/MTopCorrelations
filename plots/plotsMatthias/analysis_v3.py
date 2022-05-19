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

# from correlator_histogram_styling_v3 import style_corr_hist

import Analysis.Tools.syncer                # Starts syncing by itself, does not need to be called in script
import numpy as np

################################################################################
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory', action='store', default='MTopCorrelations_v2')
argParser.add_argument('--selection',      action='store', default='nAK81p-AK8pt')
argParser.add_argument('--era',            action='store', type=str, default="UL2018")
argParser.add_argument('--small',          action='store_true', help='Run only on a small subset of the data?', )
argParser.add_argument('--nJobs',          action='store',         nargs='?',  type=int, default=1)
argParser.add_argument('--job',            action='store',                     type=int, default=0)
argParser.add_argument('--max_cons',       action='store', type=int, default=50)
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

mc_ = [UL2018.TTbar_1]
mc_dev_ = [UL2018.TTbar_2, UL2018.TTbar_3]
if args.nJobs == 1:
    mc = mc_
    mc_dev = mc_dev_
else:
    mc = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_]
    mc_dev = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_dev_]
mc_all = mc + mc_dev
# As the #SPLIT100 function changes sample.name, only the first 11 characters of sample.name can be used!
if args.small:
    for sample in mc_all:
        sample.reduceFiles(to=1)
lumi_scale = 60

################################################################################
# Correlator Hist
hist = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist_unweighted = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist1 = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist2 = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist3 = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist4 = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist5 = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist6 = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist_top_pt = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt1 = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt2 = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt3 = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt4 = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt5 = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt6 = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_numb_all_triplets = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 100000) for sample in mc_all}
hist_numb_triplets = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 1000) for sample in mc_all}
hist_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist_unweighted_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist1_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist2_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist3_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist4_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist5_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist6_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 3) for sample in mc_all}
hist_top_pt_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt1_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt2_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt3_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt4_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt5_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_top_pt6_w = {sample.name[:11]: ROOT.TH1F("Correlator", "p_{T,top}", 40, 400, 700) for sample in mc_all}
hist_numb_all_triplets_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 100000) for sample in mc_all}
hist_numb_triplets_w = {sample.name[:11]: ROOT.TH1F("Correlator", "3 #zeta", 40, 0, 1000) for sample in mc_all}
max_numb_of_particles = args.max_cons

calc_equilateral_triplets = False
calc_isosceles_triplets = True

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
            part.SetPtEtaPhiM(event.GenJetAK8_cons_pt[i], event.GenJetAK8_cons_eta[i], event.GenJetAK8_cons_phi[i], event.GenJetAK8_cons_mass[i])
            constituents.append(part)
    def getPt(part):
        return part.Pt()
    constituents.sort(key=getPt, reverse=True)
    return constituents


def gen_tops(event, sample):
    global max_numb_of_particles, calc_equilateral_triplets, calc_isosceles_triplets

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

    event.top_had_pt = top_had.Pt()

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

    if event.merged_q_aq_b and event.nearest_jet_pt > 400:
        event.all_merged_and_high_pt = True
    else:
        event.all_merged_and_high_pt = False

    if event.all_merged_and_high_pt:
        event.matched_jet_nCons = 0
        event.matched_jet_Cons_pt = []
        for k in range(event.nGenJetAK8_cons):
            if event.GenJetAK8_cons_jetIndex[k] == nearest_jet_idx_t:
                event.matched_jet_nCons += 1
                event.matched_jet_Cons_pt.append(event.GenJetAK8_cons_pt[k])

        event.matched_jet_Cons_pt_ratio = [elem/event.nearest_jet_pt for elem in event.matched_jet_Cons_pt]
    else:
        event.matched_jet_nCons = float('nan')
        event.matched_jet_Cons_pt = float('nan')
        event.matched_jet_Cons_pt_ratio = float('nan')

    if args.nJobs > 1:
        if event.all_merged_and_high_pt:
            matched_jet_cons = getJetConstituents(event=event, idx=nearest_jet_idx_t)

            # max_numb_of_particles = 50
            if len(matched_jet_cons) < max_numb_of_particles:
                numb_of_particles = len(matched_jet_cons)
            else:
                numb_of_particles = max_numb_of_particles

            if calc_equilateral_triplets:
                numb_of_triplets = 0
                numb_of_all_triplets = 0

                delta_delta = 3.5/3. * (170./event.nearest_jet_pt)**2
                triplet = [ROOT.TLorentzVector()]*3

                for i in range(numb_of_particles):
                    for j in range(i+1, numb_of_particles):
                        for k in range(j+1, numb_of_particles):

                            numb_of_all_triplets += 1
                            triplet[0] = matched_jet_cons[i]
                            triplet[1] = matched_jet_cons[j]
                            triplet[2] = matched_jet_cons[k]
                            delta_0_1 = triplet[0].DeltaR(triplet[1])
                            delta_0_2 = triplet[0].DeltaR(triplet[2])
                            delta_1_2 = triplet[1].DeltaR(triplet[2])
                            if abs(delta_0_1-delta_0_2) < delta_delta:
                                if abs(delta_1_2-delta_0_1) < delta_delta:
                                    if abs(delta_0_2-delta_1_2) < delta_delta:
                                        numb_of_triplets += 1
                                        w = (triplet[0].Pt() * triplet[1].Pt() * triplet[2].Pt())
                                        w = w**2
                                        w = w / (event.nearest_jet_pt**(3*2))

                                        zeta = (delta_0_1 + delta_0_2 + delta_1_2) / 3

                                        hist[sample.name[:11]].Fill(zeta*3, w)
                                        hist_unweighted[sample.name[:11]].Fill(zeta*3, 1)

                                        if 400 < event.nearest_jet_pt < 450:
                                            hist1[sample.name[:11]].Fill(zeta*3, w)
                                        elif 450 < event.nearest_jet_pt < 500:
                                            hist2[sample.name[:11]].Fill(zeta*3, w)
                                        elif 500 < event.nearest_jet_pt < 550:
                                            hist3[sample.name[:11]].Fill(zeta*3, w)
                                        elif 550 < event.nearest_jet_pt < 600:
                                            hist4[sample.name[:11]].Fill(zeta*3, w)
                                        elif 600 < event.nearest_jet_pt < 650:
                                            hist5[sample.name[:11]].Fill(zeta*3, w)
                                        elif 650 < event.nearest_jet_pt < 700:
                                            hist6[sample.name[:11]].Fill(zeta*3, w)

                hist_numb_all_triplets[sample.name[:11]].Fill(numb_of_all_triplets)
                hist_numb_triplets[sample.name[:11]].Fill(numb_of_triplets)

                hist_top_pt[sample.name[:11]].Fill(event.top_had_pt)

                if 400 < event.nearest_jet_pt < 450:
                    hist_top_pt1[sample.name[:11]].Fill(event.top_had_pt)
                elif 450 < event.nearest_jet_pt < 500:
                    hist_top_pt2[sample.name[:11]].Fill(event.top_had_pt)
                elif 500 < event.nearest_jet_pt < 550:
                    hist_top_pt3[sample.name[:11]].Fill(event.top_had_pt)
                elif 550 < event.nearest_jet_pt < 600:
                    hist_top_pt4[sample.name[:11]].Fill(event.top_had_pt)
                elif 600 < event.nearest_jet_pt < 650:
                    hist_top_pt5[sample.name[:11]].Fill(event.top_had_pt)
                elif 650 < event.nearest_jet_pt < 700:
                    hist_top_pt6[sample.name[:11]].Fill(event.top_had_pt)

            if calc_isosceles_triplets:
                numb_of_triplets = 0
                numb_of_all_triplets = 0

                delta_delta = 3.5/3.*(170./event.nearest_jet_pt)**2
                short_side = 0.1
                triplet = [ROOT.TLorentzVector()]*3

                for i in range(numb_of_particles):
                    for j in range(i + 1, numb_of_particles):
                        for k in range(j + 1, numb_of_particles):

                            numb_of_all_triplets += 1
                            triplet[0] = matched_jet_cons[i]
                            triplet[1] = matched_jet_cons[j]
                            triplet[2] = matched_jet_cons[k]
                            delta_0_1 = triplet[0].DeltaR(triplet[1])
                            delta_0_2 = triplet[0].DeltaR(triplet[2])
                            delta_1_2 = triplet[1].DeltaR(triplet[2])
                            if abs(delta_0_1-delta_0_2) < delta_delta and delta_1_2 < short_side or \
                                    abs(delta_0_2-delta_1_2) < delta_delta and delta_0_1 < short_side or \
                                    abs(delta_1_2-delta_0_1) < delta_delta and delta_0_2 < short_side:
                                numb_of_triplets += 1
                                w = (triplet[0].Pt() * triplet[1].Pt() * triplet[2].Pt())
                                w = w**2
                                w = w / (event.nearest_jet_pt**(3*2))

                                zeta = (delta_0_1 + delta_0_2 + delta_1_2)/3

                                hist_w[sample.name[:11]].Fill(zeta*3, w)
                                hist_unweighted_w[sample.name[:11]].Fill(zeta*3, 1)

                                if 400 < event.nearest_jet_pt < 450:
                                    hist1_w[sample.name[:11]].Fill(zeta*3, w)
                                elif 450 < event.nearest_jet_pt < 500:
                                    hist2_w[sample.name[:11]].Fill(zeta*3, w)
                                elif 500 < event.nearest_jet_pt < 550:
                                    hist3_w[sample.name[:11]].Fill(zeta*3, w)
                                elif 550 < event.nearest_jet_pt < 600:
                                    hist4_w[sample.name[:11]].Fill(zeta*3, w)
                                elif 600 < event.nearest_jet_pt < 650:
                                    hist5_w[sample.name[:11]].Fill(zeta*3, w)
                                elif 650 < event.nearest_jet_pt < 700:
                                    hist6_w[sample.name[:11]].Fill(zeta*3, w)

                hist_numb_all_triplets_w[sample.name[:11]].Fill(numb_of_all_triplets)
                hist_numb_triplets_w[sample.name[:11]].Fill(numb_of_triplets)

                hist_top_pt_w[sample.name[:11]].Fill(event.top_had_pt)

                if 400 < event.nearest_jet_pt < 450:
                    hist_top_pt1_w[sample.name[:11]].Fill(event.top_had_pt)
                elif 450 < event.nearest_jet_pt < 500:
                    hist_top_pt2_w[sample.name[:11]].Fill(event.top_had_pt)
                elif 500 < event.nearest_jet_pt < 550:
                    hist_top_pt3_w[sample.name[:11]].Fill(event.top_had_pt)
                elif 550 < event.nearest_jet_pt < 600:
                    hist_top_pt4_w[sample.name[:11]].Fill(event.top_had_pt)
                elif 600 < event.nearest_jet_pt < 650:
                    hist_top_pt5_w[sample.name[:11]].Fill(event.top_had_pt)
                elif 650 < event.nearest_jet_pt < 700:
                    hist_top_pt6_w[sample.name[:11]].Fill(event.top_had_pt)


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

for sample in mc_dev:
    sample.style = styles.lineStyle(sample.color)

for sample in mc_all:
    sample.weight = lambda event, sample: 1.

stack = Stack(mc, [], mc_dev)

# Use some defaults
Plot.setDefaults(stack = stack, weight = staticmethod(weight_), selectionString = cutInterpreter.cutString(args.selection))

################################################################################
# Now define the plots

plots = []

plots.append(Plot(
    name="Dummy_Plot",
    texX='Dummy plot to make plotting.fill() work', texY='Number of Events',
    attribute=lambda event, sample: 1.,
    binning=[3, -0.5, 10.5],
))

if args.nJobs == 1:
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
        attribute = lambda event, sample: event.delta_min if event.all_merged_and_high_pt else float('nan'),
        binning=[25, 0., 7.],
    ))

    plots.append(Plot(
        name = "cut_jet_mass_2",
        texX = 'Mass of jet next to hadronic top [GeV]', texY = 'Number of Events',
        attribute = lambda event, sample: event.nearest_jet_mass_t if event.all_merged_and_high_pt else float('nan'),
        binning=[25, 0., 500.],
    ))

    plots.append(Plot(
        name = "cut_jet_pt_2",
        texX = 'P_{t} of jet next to hadronic top [GeV]', texY = 'Number of Events',
        attribute = lambda event, sample: event.nearest_jet_pt if event.all_merged_and_high_pt else float('nan'),
        binning=[25, 0., 500.],
    ))

    plots.append(Plot(
        name = "cut_jet_eta_2",
        texX = '#eta of jet next to hadronic top [GeV]', texY = 'Number of Events',
        attribute = lambda event, sample: event.nearest_jet_eta if event.all_merged_and_high_pt else float('nan'),
        binning=[25, 3.5, 3.5],
    ))

    plots.append(Plot(
        name = "cut_jet_phi_2",
        texX = '#phi of jet next to hadronic top [GeV]', texY = 'Number of Events',
        attribute = lambda event, sample: event.nearest_jet_phi if event.all_merged_and_high_pt else float('nan'),
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

    plots.append(Plot(
        name = "ratio_top_pt_jet_pt",
        texX = 'Ratio of p_{T,top} and p_{T,jet}', texY = 'Number of Events',
        attribute = lambda event, sample: event.top_had_pt/event.nearest_jet_pt,
        binning=[25, 0.8, 1.2],
    ))


plotting.fill(plots, read_variables = read_variables, sequence = sequence)

if args.nJobs == 1:
    drawPlots(plots)

# print hist.Integral()

f = ROOT.TFile('correlator_part_{:}_{:}.root'.format(max_numb_of_particles, args.job), 'RECREATE')
f.cd()
if calc_equilateral_triplets:
    for sample_name, single_hist in hist.items(): single_hist.Write('correlator_hist'+sample_name)
    for sample_name, single_hist in hist_unweighted.items(): single_hist.Write('correlator_hist_unweighted'+sample_name)
    for sample_name, single_hist in hist1.items(): single_hist.Write('correlator_hist_400_450'+sample_name)
    for sample_name, single_hist in hist2.items(): single_hist.Write('correlator_hist_450_500'+sample_name)
    for sample_name, single_hist in hist3.items(): single_hist.Write('correlator_hist_500_550'+sample_name)
    for sample_name, single_hist in hist4.items(): single_hist.Write('correlator_hist_550_600'+sample_name)
    for sample_name, single_hist in hist5.items(): single_hist.Write('correlator_hist_600_650'+sample_name)
    for sample_name, single_hist in hist6.items(): single_hist.Write('correlator_hist_650_700'+sample_name)
    for sample_name, single_hist in hist_top_pt.items(): single_hist.Write('top_pt_hist'+sample_name)
    for sample_name, single_hist in hist_top_pt1.items(): single_hist.Write('top_pt_hist_400_450'+sample_name)
    for sample_name, single_hist in hist_top_pt2.items(): single_hist.Write('top_pt_hist_450_500'+sample_name)
    for sample_name, single_hist in hist_top_pt3.items(): single_hist.Write('top_pt_hist_500_550'+sample_name)
    for sample_name, single_hist in hist_top_pt4.items(): single_hist.Write('top_pt_hist_550_600'+sample_name)
    for sample_name, single_hist in hist_top_pt5.items(): single_hist.Write('top_pt_hist_600_650'+sample_name)
    for sample_name, single_hist in hist_top_pt6.items(): single_hist.Write('top_pt_hist_650_700'+sample_name)
    for sample_name, single_hist in hist_numb_all_triplets.items(): single_hist.Write('number_of_all_triplets'+sample_name)
    for sample_name, single_hist in hist_numb_triplets.items(): single_hist.Write('number_of_selected_triplets'+sample_name)

if calc_isosceles_triplets:
    for sample_name, single_hist in hist_w.items(): single_hist.Write('correlator_hist_w'+sample_name)
    for sample_name, single_hist in hist_unweighted_w.items(): single_hist.Write('correlator_hist_unweighted_w'+sample_name)
    for sample_name, single_hist in hist1_w.items(): single_hist.Write('correlator_hist_400_450_w'+sample_name)
    for sample_name, single_hist in hist2_w.items(): single_hist.Write('correlator_hist_450_500_w'+sample_name)
    for sample_name, single_hist in hist3_w.items(): single_hist.Write('correlator_hist_500_550_w'+sample_name)
    for sample_name, single_hist in hist4_w.items(): single_hist.Write('correlator_hist_550_600_w'+sample_name)
    for sample_name, single_hist in hist5_w.items(): single_hist.Write('correlator_hist_600_650_w'+sample_name)
    for sample_name, single_hist in hist6_w.items(): single_hist.Write('correlator_hist_650_700_w'+sample_name)
    for sample_name, single_hist in hist_top_pt_w.items(): single_hist.Write('top_pt_hist_w'+sample_name)
    for sample_name, single_hist in hist_top_pt1_w.items(): single_hist.Write('top_pt_hist_400_450_w'+sample_name)
    for sample_name, single_hist in hist_top_pt2_w.items(): single_hist.Write('top_pt_hist_450_500_w'+sample_name)
    for sample_name, single_hist in hist_top_pt3_w.items(): single_hist.Write('top_pt_hist_500_550_w'+sample_name)
    for sample_name, single_hist in hist_top_pt4_w.items(): single_hist.Write('top_pt_hist_550_600_w'+sample_name)
    for sample_name, single_hist in hist_top_pt5_w.items(): single_hist.Write('top_pt_hist_600_650_w'+sample_name)
    for sample_name, single_hist in hist_top_pt6_w.items(): single_hist.Write('top_pt_hist_650_700_w'+sample_name)
    for sample_name, single_hist in hist_numb_all_triplets_w.items(): single_hist.Write('number_of_all_triplets_w'+sample_name)
    for sample_name, single_hist in hist_numb_triplets_w.items(): single_hist.Write('number_of_selected_triplets_w'+sample_name)

f.Close()

# for pt_range in ['', '_unweighted', '_400_450', '_450_500', '_500_550', '_550_600', '_600_650', '_650_700']:
#     style_corr_hist(filename_root='correlator_'+str(args.job)+'.root',
#                     hist_name='correlator_hist'+pt_range,
#                     filename_graphic='/correlator_hist'+pt_range+'.png')

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

# Root files zusammen addieren mit Konsolen-Befehl:
# hadd correlator-joined.root correlator_1.root correlator_2.root correlator_3.root
# oder
# hadd correlator-joined.root correlator_*.root

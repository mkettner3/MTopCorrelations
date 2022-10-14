# -*- coding: utf-8 -*-

"""
Script to calculate triplets and generate histograms event-wise. The triplet data is only stored in RAM for one event.
"""

import time
from math import isnan
import MTopCorrelations.samples.nanoTuples_UL2018_nanoAOD as UL2018
from triplet_maker import make_triplets_and_cut
from find_hadronic_jet import find_hadronic_jet
from jet_constituents import get_jet_constituents
import ROOT
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse


def calc_triplets_and_hist(samples, pt_jet_ranges, max_delta_zeta=float('nan'), delta_legs=float('nan'),
                           shortest_side=0.1, nbins=50, hist_range=(0, 3), pt_variations=None):
    global count, number_events

    read_variables = [
        "Generator_weight/D",
        "nGenPart/I",
        "GenPart[pt/F,eta/F,phi/F,m/F,pdgId/I,mompdgId/I,grmompdgId/I]",
        "nGenJetAK8/I",
        "GenJetAK8[pt/F,eta/F,phi/F,mass/F]",
        "nGenJetAK8_cons/I",
        VectorTreeVariable.fromString("GenJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000),
        "nPFJetAK8/I",
        "PFJetAK8[pt/F,eta/F,phi/F,mass/F]",
        "nPFJetAK8_cons/I",
        VectorTreeVariable.fromString("PFJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000)]

    hists = [[[[ROOT.TH1F("Correlator", "3 #zeta", nbins, hist_range[0], hist_range[1]) for _ in range(2)] for _ in range(len(pt_jet_ranges))] for _ in range(len(samples))] for _ in range(2)]
    hists_w = [[[[ROOT.TH1F("Correlator", "3 #zeta", nbins, hist_range[0], hist_range[1]) for _ in range(2)] for _ in range(len(pt_jet_ranges))] for _ in range(len(samples))] for _ in range(2)]
    hists_varied = [[[[ROOT.TH1F("Correlator", "3 #zeta", nbins, hist_range[0], hist_range[1]) for _ in range(2)] for _ in range(len(pt_jet_ranges))] for _ in range(len(samples))] for _ in range(2)]
    hists_jet_pt = [ROOT.TH1F("Hadronic Top-Jet-p_{t}", "Jet-p_{t}", nbins, 380, 730) for _ in range(2)]
    hists_jet_mass = [ROOT.TH1F("Hadronic Top-Jet-mass", "Jet-mass", nbins, 75, 300) for _ in range(2)]
    hists_event_weight = [[ROOT.TH1F("Event weights", "event-weight", nbins, -1, 1) for _ in range(len(samples))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample in enumerate(samples):
            r = sample.treeReader(variables=read_variables, selectionString="Sum$({:}JetAK8_pt>400)>=1".format(level))
            r.start()
            while r.run():                                                              # Event-Loop
                event_weight = 60 * 831.762 * 3*0.108 * (1-3*0.108)*2 * 1000 / number_events[h] * r.event.Generator_weight

                hadronic_jet_idx, hadronic_jet_pt, hadronic_jet_mass = find_hadronic_jet(r.event, level=level,
                                                                                         merge_tolerance=0.8,
                                                                                         jet_pt_min=400)

                if hadronic_jet_idx is not None:
                    jet_constituents = get_jet_constituents(event=r.event, level=level,
                                                            index=hadronic_jet_idx, particle_types='charged_par',
                                                            max_numb_of_cons=50)

                    if len(jet_constituents) > 0:
                        triplets, triplets_w = construct_triplets(hadronic_jet_pt, jet_constituents,
                                                                  max_delta_zeta, delta_legs, shortest_side)

                        for k, jet_range in enumerate(pt_jet_ranges):
                            if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                for three_zeta, weight in zip(triplets[0], triplets[1]):
                                    hists[g][h][k][0].Fill(three_zeta, weight*event_weight)
                                    hists[g][h][k][1].Fill(three_zeta, event_weight)
                                for three_zeta, weight in zip(triplets_w[0], triplets_w[1]):
                                    hists_w[g][h][k][0].Fill(three_zeta, weight*event_weight)
                                    hists_w[g][h][k][1].Fill(three_zeta, event_weight)
                                break

                        hists_jet_pt[g].Fill(hadronic_jet_pt, event_weight)
                        hists_jet_mass[g].Fill(hadronic_jet_mass, event_weight)
                        hists_event_weight[g][h].Fill(event_weight)
                        count += 1

                        for v, var_fac in enumerate(pt_variations):
                            hadronic_jet_pt, jet_constituents = vary_pt(jet_pt=hadronic_jet_pt,
                                                                        constituents=jet_constituents, factor=var_fac)

                            triplets, triplets_w = construct_triplets(hadronic_jet_pt, jet_constituents,
                                                                      max_delta_zeta, delta_legs, shortest_side)

                            for k, jet_range in enumerate(pt_jet_ranges):
                                if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                    for three_zeta, weight in zip(triplets[0], triplets[1]):
                                        hists_varied[g][h][k][v].Fill(three_zeta, weight*event_weight)
                                    break

    return hists, hists_w, hists_jet_pt, hists_jet_mass, hists_varied, hists_event_weight


def construct_triplets(hadronic_jet_pt, jet_constituents, max_delta_zeta, delta_legs, shortest_side):
    # type: (float, list, float, float, float) -> tuple

    if isnan(max_delta_zeta):
        max_delta_zeta_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
    else:
        max_delta_zeta_calc = max_delta_zeta

    if isnan(delta_legs):
        delta_legs_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
    else:
        delta_legs_calc = delta_legs

    (triplets,
     triplets_w) = make_triplets_and_cut(jet_pt=hadronic_jet_pt, particle_vectors=jet_constituents,
                                         max_delta_zeta=max_delta_zeta_calc,
                                         delta_legs=delta_legs_calc, shortest_side=shortest_side)

    return triplets, triplets_w


def vary_pt(jet_pt, constituents, factor):
    # type: (float, list, float) -> tuple

    jet_pt = jet_pt * factor
    # for i in range(len(constituents)):
    #     constituents[i] = constituents[i] * factor

    return jet_pt, constituents


def save_root_hists(hists_top, hists_w, hists_jet_pt, hists_jet_mass, hists_varied, pt_variations, hists_ev_weight,
                    sample_names, pt_jet_ranges, filename):
    # type: (list, list, list, list, list, tuple, list, list, list, str) -> None

    r_file = ROOT.TFile(filename, 'RECREATE')
    r_file.cd()
    for p_type in ['Top-Quark', 'W-Boson', 'Others']:
        for lev in ['Gen-Level', 'PF-Level']:
            if p_type in ['Top-Quark', 'W-Boson']:
                for weight in ['weighted', 'absolute']:
                    r_file.mkdir(p_type+'/'+lev+'/'+weight)
            else:
                r_file.mkdir(p_type+'/'+lev)

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            for k, pt_jet_range in enumerate(pt_jet_ranges):
                for f, (weighted, weight_dir) in enumerate(zip(['', '_abscou'], ['weighted', 'absolute'])):
                    r_file.cd('/Top-Quark/'+level+'-Level/'+weight_dir)
                    hists_top[g][h][k][f].Write('correlator_hist_{:}_{:}_{:}_{:}{:}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1], weighted))
                    r_file.cd('/W-Boson/'+level+'-Level/'+weight_dir)
                    hists_w[g][h][k][f].Write('correlator_hist_W_{:}_{:}_{:}_{:}{:}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1], weighted))
                for v, var_fac in enumerate(pt_variations):
                    r_file.cd('/Top-Quark/'+level+'-Level/weighted')
                    hists_varied[g][h][k][v].Write('correlator_hist_varied_{:.2f}_{:}_{:}_{:}_{:}'.format(var_fac, level, sample_name, pt_jet_range[0], pt_jet_range[1]))
        r_file.cd('/Others/'+level+'-Level')
        hists_jet_pt[g].Write('hadronic_top_jet_pt_hist_{:}'.format(level))
        hists_jet_mass[g].Write('hadronic_top_jet_mass_hist_{:}'.format(level))
        for h, sample_name in enumerate(sample_names):
            hists_ev_weight[g][h].Write('event_weights_{:}_{:}'.format(level, sample_name))

    r_file.Close()


mc_ttbar = [UL2018.TTbar_1, UL2018.TTbar_2, UL2018.TTbar_3, UL2018.TTbar_4, UL2018.TTbar_5]
mc_qcd = [UL2018.QCD_0]
number_events_ttbar = [61694526814.8, 55968875821.1, 139438020309, 51889503097.9, 42054325908.3]
number_events_qcd = [13849346.0]

pt_jet_lowest = 400
pt_jet_highest = 700
pt_jet_step = 50
pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                    range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    argParser.add_argument('--sample_type', action='store', type=str, default='TTbar')  # Can either be 'TTbar' or 'QCD'
    args = argParser.parse_args()

    mc = mc_ttbar if args.sample_type != 'QCD' else mc_qcd
    number_events = number_events_ttbar if args.sample_type != 'QCD' else number_events_qcd
    samples = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc]

    start = time.time()
    count = 0
    (hists, hists_w,
     hists_jet_pt, hists_jet_mass,
     hists_varied,
     hists_event_weight) = calc_triplets_and_hist(samples=samples, pt_jet_ranges=pt_jet_ranges,
                                                  max_delta_zeta=float('nan'), delta_legs=float('nan'), shortest_side=0.05,
                                                  pt_variations=(1.02, 0.98))
    save_root_hists(hists_top=hists, hists_w=hists_w, hists_jet_pt=hists_jet_pt, hists_jet_mass=hists_jet_mass,
                    hists_varied=hists_varied, pt_variations=(1.02, 0.98), hists_ev_weight=hists_event_weight,
                    sample_names=[sample.name[:11] for sample in samples], pt_jet_ranges=pt_jet_ranges,
                    filename='histogram_files/correlator_hist_trip_13_pp_{:03}_{:}.root'.format(args.job, args.sample_type))
    end = time.time()

    print('Executing calc_triplet_and_hist.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))
    print('Number of considered events in all samples: {:}'.format(count))

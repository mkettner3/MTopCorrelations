# -*- coding: utf-8 -*-

"""
Script to calculate triplets and generate histograms event-wise. The triplet data is only stored in RAM for one event.
"""

import time
from math import isnan
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
from MTopCorrelations.Tools.triplet_maker import make_triplets_and_cut
from MTopCorrelations.Tools.jet_constituents import get_jet_constituents
from calc_triplet_data import find_hadronic_jet
import ROOT
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse


def calc_triplets_and_hist(samples, pt_jet_ranges, max_delta_zeta=float('nan'), delta_legs=float('nan'), shortest_side=0.1,
                           nbins=50, hist_range=(0, 3)):
    read_variables = [
        "Generator_weight/F",
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
    hists_jet_pt = [ROOT.TH1F("Hadronic Top-Jet-p_{t}", "Jet-p_{t}", nbins, 380, 730) for _ in range(2)]
    hists_jet_mass = [ROOT.TH1F("Hadronic Top-Jet-mass", "Jet-mass", nbins, 0, 450) for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample in enumerate(samples):
            r = sample.treeReader(variables=read_variables, selectionString="Sum$({:}JetAK8_pt>400)>=1".format(level))

            global count
            # count_int = 0
            r.start()
            while r.run():                                                              # Event-Loop
                event_weight = 60 * 831.762 * 3*0.108 * (1-3*0.108)*2 * 1000 / number_events[h] * r.event.Generator_weight

                # if count_int < 10:
                #     print('Event-Weight ({:}): {:}'.format(sample.name, event_weight))

                hadronic_jet_idx, hadronic_jet_pt, hadronic_jet_mass = find_hadronic_jet(r.event, level=level,
                                                                                         merge_tolerance=0.8,
                                                                                         jet_pt_min=400)

                if hadronic_jet_idx is not None:
                    jet_constituents = get_jet_constituents(event=r.event, level=level,
                                                            index=hadronic_jet_idx, max_numb_of_cons=50)

                    if len(jet_constituents) > 0:
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

                        for k, jet_range in enumerate(pt_jet_ranges):
                            if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                for three_zeta, weight in zip(triplets[0], triplets[1]):
                                    hists[g][h][k][0].Fill(three_zeta, weight*event_weight)
                                    hists[g][h][k][1].Fill(three_zeta, event_weight)
                                for three_zeta, weight in zip(triplets_w[0], triplets_w[1]):
                                    hists_w[g][h][k][0].Fill(three_zeta, weight*event_weight)
                                    hists_w[g][h][k][1].Fill(three_zeta, event_weight)
                                break

                        # if count_int < 10:
                        #     print('Corr-Weight ({:}): {:}'.format(sample.name, weight))

                        hists_jet_pt[g].Fill(hadronic_jet_pt, event_weight)
                        hists_jet_mass[g].Fill(hadronic_jet_mass, event_weight)
                        count += 1
                        # count_int += 1

    return hists, hists_w, hists_jet_pt, hists_jet_mass


def save_root_hists(hists_top, hists_w, hists_jet_pt, hists_jet_mass, sample_names, pt_jet_ranges, filename):
    # type: (list, list, list, list, list, list, str) -> None

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
        r_file.cd('/Others/'+level+'-Level')
        hists_jet_pt[g].Write('hadronic_top_jet_pt_hist_{:}'.format(level))
        hists_jet_mass[g].Write('hadronic_top_jet_mass_hist_{:}'.format(level))

    r_file.Close()


mc = [UL2018.TTbar_1, UL2018.TTbar_2, UL2018.TTbar_3]
number_events = [14815731313, 17230816904, 12779934801]

pt_jet_lowest = 400
pt_jet_highest = 700
pt_jet_step = 50
pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                    range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    args = argParser.parse_args()

    samples = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc]

    start = time.time()
    count = 0
    hists, hists_w, hists_jet_pt, hists_jet_mass = calc_triplets_and_hist(samples=samples, pt_jet_ranges=pt_jet_ranges,
                                                                          max_delta_zeta=float('nan'),
                                                                          delta_legs=float('nan'), shortest_side=0.05)
    save_root_hists(hists_top=hists, hists_w=hists_w, hists_jet_pt=hists_jet_pt, hists_jet_mass=hists_jet_mass,
                    sample_names=[sample.name[:11] for sample in samples],
                    pt_jet_ranges=pt_jet_ranges,
                    filename='histogram_files/correlator_hist_trip_pp_{:02}.root'.format(args.job))
    end = time.time()

    print('Executing calc_triplet_and_hist.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))
    print('Number of considered events in all samples: {:}'.format(count))

# -*- coding: utf-8 -*-

"""
Script to calculate triplets and generate histograms event-wise. The triplet data is only stored in RAM for one event.
"""

import numpy as np
import time
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
from MTopCorrelations.Tools.triplet_maker import make_triplets_and_cut
from MTopCorrelations.Tools.jet_constituents import get_jet_constituents
from calc_triplet_data import find_hadronic_jet
from gen_hist_from_triplets import store_np_hist_in_root
import ROOT
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse


def calc_triplets_and_hist(samples, pt_jet_ranges, max_delta_zeta=None, delta_legs=None, shortest_side=None,
                           nbins=50, hist_range=(0, 3)):
    read_variables = [
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

    np_hists = [[[np.full(nbins, 0, dtype=np.float32) for _ in range(len(pt_jet_ranges))] for _ in range(len(samples))] for _ in range(2)]
    # hists = [[ROOT.TH1F("Correlator", "3 #zeta", nbins, hist_range[0], hist_range[1]) for _ in range(len(pt_jet_ranges))] for _ in range(len(samples))]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample in enumerate(samples):
            r = sample.treeReader(variables=read_variables, selectionString="Sum$({:}JetAK8_pt>400)>=1".format(level))

            global count
            r.start()
            while r.run():                                                              # Event-Loop
                hadronic_jet_idx, hadronic_jet_pt = find_hadronic_jet(r.event, level=level,
                                                                      merge_tolerance=0.8, jet_pt_min=400)

                if hadronic_jet_idx is not None:
                    jet_constituents = get_jet_constituents(event=r.event, level=level,
                                                            index=hadronic_jet_idx, max_numb_of_cons=50)

                    if len(jet_constituents) > 0:
                        if max_delta_zeta is None:
                            max_delta_zeta_calc = max_delta_zeta
                        elif np.isnan(max_delta_zeta):
                            max_delta_zeta_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
                        else:
                            max_delta_zeta_calc = max_delta_zeta

                        if delta_legs is None:
                            delta_legs_calc = delta_legs
                        elif np.isnan(delta_legs):
                            delta_legs_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
                        else:
                            delta_legs_calc = delta_legs

                        triplets = make_triplets_and_cut(jet_pt=hadronic_jet_pt, particle_vectors=jet_constituents,
                                                         max_delta_zeta=max_delta_zeta_calc,
                                                         delta_legs=delta_legs_calc, shortest_side=shortest_side)

                        k = None
                        for i, jet_range in enumerate(pt_jet_ranges):
                            if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                # for three_zeta, weight in zip(triplets[0], triplets[1]):
                                #     hists[h][k].Fill(three_zeta, weight)
                                k = i
                                break

                        if k is not None:
                            np_hists[g][h][k] = np_hists[g][h][k] + np.histogram(triplets[0], weights=triplets[1],
                                                                                 bins=nbins, range=hist_range,
                                                                                 density=False)[0]
                            # ToDo: Test if it is faster to directly fill the root histogram
                        count += 1

    np_hist_bins = np.histogram(triplets[0], weights=triplets[1],
                                bins=nbins, range=hist_range, density=False)[1]

    return np_hists, np_hist_bins


def save_root_hists(root_hists, sample_names, pt_jet_ranges, filename):
    # type: (list, list, list, str) -> None

    f = ROOT.TFile(filename, 'RECREATE')
    f.cd()

    for h, sample_name in enumerate(sample_names):
        for k, pt_jet_range in enumerate(pt_jet_ranges):
            root_hists[h][k].Write('correlator_hist_{:}_{:}_{:}'.format(sample_name, pt_jet_range[0], pt_jet_range[1]))

    f.Close()


mc = [UL2018.TTbar_1, UL2018.TTbar_2, UL2018.TTbar_3]

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
    numpy_hist = calc_triplets_and_hist(samples=samples, pt_jet_ranges=pt_jet_ranges, max_delta_zeta=np.nan)
                                        # delta_legs=np.nan, shortest_side=0.1)
    # save_root_hists(root_hists=hists, sample_names=[sample.name[:11] for sample in samples],
    #                 pt_jet_ranges=pt_jet_ranges,
    #                 filename='histogram_files/correlator_hist_trip_test_new_{:02}.root'.format(args.job))
    store_np_hist_in_root(numpy_hist=numpy_hist, sample_names=[sample.name[:11] for sample in samples],
                          pt_jet_ranges=pt_jet_ranges,
                          filename='histogram_files/correlator_hist_trip_pp_{:02}.root'.format(args.job))
    end = time.time()

    print('Executing calc_triplet_and_hist.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))
    print('Number of considered events in all samples: {:}'.format(count))

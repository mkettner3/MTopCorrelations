# -*- coding: utf-8 -*-

"""
Script to calculate triplets and generate histograms event-wise. The triplet data is only stored in RAM for one event.
"""

import numpy as np
import time
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
from MTopCorrelations.Tools.triplet_maker import make_triplets
from MTopCorrelations.Tools.jet_constituents import get_jet_constituents
from calc_triplet_data import find_hadronic_jet
from gen_hist_from_triplets import store_np_hist_in_root
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse


def calc_triplets_and_hist(samples, pt_jet_ranges, max_delta_zeta=None, delta_legs=None, shortest_side=None,
                           nbins=40, hist_range=(0, 3)):
    read_variables = [
        "nGenPart/I",
        "GenPart[pt/F,eta/F,phi/F,m/F,pdgId/I,mompdgId/I,grmompdgId/I]",
        "nGenJetAK8/I",
        "GenJetAK8[pt/F,eta/F,phi/F,mass/F]",
        "nGenJetAK8_cons/I",
        VectorTreeVariable.fromString("GenJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000)]

    np_hists = [[np.full(nbins, 0, dtype=np.float32)]*len(pt_jet_ranges)]*len(samples)

    for h, sample in enumerate(samples):
        r = sample.treeReader(variables=read_variables, selectionString="Sum$(GenJetAK8_pt>400)>=1")

        global count
        r.start()
        while r.run():                                                              # Event-Loop
            # count += 1
            # if count % 100 == 0:
            #     print('Event {} is calculated.'.format(count))
            hadronic_jet_idx, hadronic_jet_pt = find_hadronic_jet(r.event, merge_tolerance=0.8, jet_pt_min=400)

            if hadronic_jet_idx is not None:
                jet_constituents = get_jet_constituents(event=r.event, index=hadronic_jet_idx, max_numb_of_cons=50)

                if len(jet_constituents) > 0:
                    if count < 10:                      # TEST-Line
                        print(hadronic_jet_idx)         # TEST-Line
                    triplets = make_triplets(jet_pt=hadronic_jet_pt, particle_vectors=jet_constituents,
                                             top_data=bool(max_delta_zeta), w_data=(bool(delta_legs) and bool(shortest_side)),
                                             pt_value=False)

                    if np.isnan(max_delta_zeta):
                        max_delta_zeta_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
                    else:
                        max_delta_zeta_calc = max_delta_zeta

                    if max_delta_zeta is not None:
                        triplets_cut = (triplets[:, 2] < max_delta_zeta_calc)
                        if delta_legs is not None and shortest_side is not None:
                            triplets_cut = triplets_cut & (triplets[:, 3] < delta_legs) & (triplets[:, 4] < shortest_side)
                    elif delta_legs is not None and shortest_side is not None:
                        triplets_cut = (triplets[:, 2] < delta_legs) & (triplets[:, 3] < shortest_side)
                    else:
                        triplets_cut = np.full(triplets.shape[0], True, dtype=bool)

                    k = None
                    for i, jet_range in enumerate(pt_jet_ranges):
                        if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                            k = i
                            break

                    if k is not None:
                        np_hists[h][k] = np_hists[h][k] + np.histogram(triplets[triplets_cut, 0],
                                                                       weights=triplets[triplets_cut, 1],
                                                                       bins=nbins, range=hist_range, density=False)[0]
                        # ToDo: Test if it is faster to directly fill the root histogram
                    count += 1

    np_hist_bins = np.histogram(triplets[triplets_cut, 0], weights=triplets[triplets_cut, 1],
                                bins=nbins, range=hist_range, density=False)[1]

    return np_hists, np_hist_bins


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    args = argParser.parse_args()

    mc_ = [UL2018.TTbar_1]            # + [UL2018.TTbar_2, UL2018.TTbar_3]
    samples = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_]

    pt_jet_lowest = 400
    pt_jet_highest = 700
    pt_jet_step = 50
    pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                        range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))

    start = time.time()
    count = 0
    np_hists, np_hist_bins = calc_triplets_and_hist(samples=samples, pt_jet_ranges=pt_jet_ranges,
                                                    max_delta_zeta=np.nan)
    store_np_hist_in_root(numpy_hist=(np_hists, np_hist_bins), sample_names=[sample.name[:11] for sample in samples],
                          pt_jet_ranges=pt_jet_ranges,
                          filename='histogram_files/correlator_hist_trip_{:02}.root'.format(args.job))
    end = time.time()

    print('Executing calc_triplet_and_hist.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))
    print('Number of considered events: {:}'.format(count))

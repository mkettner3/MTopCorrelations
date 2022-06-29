# -*- coding: utf-8 -*-

"""
Script to generate histograms from a numpy array of triplet parameters.
"""

import numpy as np
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
import ROOT
import argparse


def read_triplets_from_hdf5_file(filename):
    # type: (str) -> np.ndarray

    import h5py

    with h5py.File(filename, 'r') as f:
        all_triplets = f['triplet-parameters'][:, :]

    return all_triplets


def calc_triplet_histogram(triplets, jet_pt_range, max_delta_zeta=None, delta_legs=None, shortest_side=None,
                           nbins=40, hist_range=(0, 3)):
    # type: (np.ndarray, tuple, float, float, float, int, tuple) -> tuple

    if max_delta_zeta is not None:
        triplets_cut = (triplets[:, 2] < max_delta_zeta) & (triplets[:, 5] > jet_pt_range[0]) & (triplets[:, 5] < jet_pt_range[1])
    elif delta_legs is not None and shortest_side is not None:
        triplets_cut = (triplets[:, 3] < delta_legs) & (triplets[:, 4] < shortest_side) & \
                       (triplets[:, 5] > jet_pt_range[0]) & (triplets[:, 5] < jet_pt_range[1])
    else:
        triplets_cut = np.full(triplets.shape[0], True, dtype=bool)

    np_hist = np.histogram(triplets[triplets_cut, 0], weights=triplets[triplets_cut, 1], bins=nbins, range=hist_range,
                           density=False)

    return np_hist


def store_np_hist_in_root(numpy_hist, sample_names, pt_jet_ranges, filename):
    # type: (tuple, list, list, str) -> None

    f = ROOT.TFile(filename, 'RECREATE')
    f.cd()

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            for j, pt_jet_range in enumerate(pt_jet_ranges):
                hist = ROOT.TH1F("Correlator", "3 #zeta", len(numpy_hist[0][g][h][j]), numpy_hist[1][0], numpy_hist[1][-1])
                for i, hist_value in enumerate(numpy_hist[0][g][h][j]):
                    hist.Fill(np.mean([numpy_hist[1][i+1], numpy_hist[1][i]]), hist_value)

                hist.Write('correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1]))
                del hist

    f.Close()


def store_np_hist_in_root_old(numpy_hist, pt_jet_ranges, filename):
    # type: (tuple, list, str) -> None

    f = ROOT.TFile(filename, 'RECREATE')
    f.cd()

    for j, pt_jet_range in enumerate(pt_jet_ranges):
        hist = ROOT.TH1F("Correlator", "3 #zeta", len(numpy_hist[0][j]), numpy_hist[1][0], numpy_hist[1][-1])
        for i, hist_value in enumerate(numpy_hist[0][j]):
            hist.Fill(np.mean([numpy_hist[1][i+1], numpy_hist[1][i]]), hist_value)

        hist.Write('correlator_hist_{:}_{:}'.format(pt_jet_range[0], pt_jet_range[1]))

    f.Close()


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    args = argParser.parse_args()

    mc_ = [UL2018.TTbar_1]            # + [UL2018.TTbar_2, UL2018.TTbar_3]
    samples = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_]

    pt_jet_range = (450, 500)

    for sample in samples:
        triplets = read_triplets_from_hdf5_file(filename='triplet_files/EWC_triplets_{:}_{:02}.h5'.format(sample.name[:11], args.job))
        delta_zeta = 3.5 / 3. * (170. / np.mean(triplets[:, 5])) ** 2
        np_hist = calc_triplet_histogram(triplets=triplets, jet_pt_range=pt_jet_range, max_delta_zeta=delta_zeta)
        store_np_hist_in_root_old(numpy_hist=np_hist, pt_jet_ranges=[pt_jet_range],
                                  filename='histogram_files/correlator_hist_trip_{:}_{:02}.root'.format(sample.name[:11],
                                                                                                        args.job))

# -*- coding: utf-8 -*-

"""
Function to normalize the covariance matrix
"""

from typing import Any
import numpy as np
import ROOT


def calc_norm_cov_matrix(filename_root_hist, hist_name):
    # type: (str, str) -> np.ndarray
    """
    Documentation
    """

    f = ROOT.TFile(filename_root_hist)
    root_hist = f.Get(hist_name)

    num_bins = root_hist.GetNbinsX()
    matrix_orig = np.zeros((num_bins, num_bins), dtype=np.float64)
    for i in range(num_bins):
        matrix_orig[i, i] = root_hist.GetBinError(i+1)**2

    matrix_norm = normalize_cov_matrix(matrix_orig=matrix_orig, root_hist=root_hist)

    return matrix_norm


def normalize_cov_matrix(matrix_orig, root_hist):
    # type: (np.ndarray, Any) -> np.ndarray

    num_bins = matrix_orig.shape[0]
    integral = root_hist.Integral()

    matrix_norm = np.zeros((num_bins, num_bins), dtype=np.float64)
    for i in range(num_bins):
        for j in range(num_bins):
            for idx in [i, j]:
                bincontent = root_hist.GetBinContent(idx)
                binwidth = root_hist.GetBinWidth(idx)
                derivation_a = - bincontent / (integral ** 2 * binwidth)
                derivation_b = (integral - bincontent) / (integral ** 2 * binwidth)
                deriv_mat = np.full((num_bins, num_bins), derivation_a, dtype=np.float64)
                if idx == i:
                    deriv_mat[i, :] = derivation_b
                elif idx == j:
                    deriv_mat[:, j] = derivation_b
                matrix_orig = deriv_mat * matrix_orig
            matrix_norm[i, j] = matrix_orig.sum()

    return matrix_norm


def normalize_cov_matrix_dennis(matrix_orig, root_hist):
    # type: (np.ndarray, Any) -> np.ndarray

    num_bins = matrix_orig.shape[0]
    integral = 10

    matrix_norm = np.zeros((num_bins, num_bins), dtype=np.float64)
    for i in range(num_bins):
        for j in range(num_bins):
            bincontent_i = 5                   # root_hist.GetBinContent(idx+1)
            bincontent_j = 5
            binwidth_i = 1                     # root_hist.GetBinWidth(idx+1)
            binwidth_j = 1
            derivation_i_a = - bincontent_i / (integral ** 2 * binwidth_i)
            derivation_i_b = (integral - bincontent_i) / (integral ** 2 * binwidth_i)
            derivation_j_a = - bincontent_j / (integral ** 2 * binwidth_j)
            derivation_j_b = (integral - bincontent_j) / (integral ** 2 * binwidth_j)
            sum_element = 0
            for h in range(num_bins):
                for k in range(num_bins):
                    derivation_i = derivation_i_b if i == k else derivation_i_a
                    derivation_j = derivation_j_b if j == h else derivation_j_a
                    sum_element += derivation_i * derivation_j * matrix_orig[h, k]
            matrix_norm[i, j] = sum_element

    return matrix_norm


def NormalizeMatrix(old_cov, hist_):

    nbins = old_cov.shape[0]
    new_cov = np.empty((nbins, nbins), dtype=np.float64)
    integral = hist_.Integral()
    for i in range(1, nbins+1):
        for j in range(1, nbins+1):
            sum_value = 0
            for k in range(1, nbins+1):
                for l in range(1, nbins+1):
                    old_entry = old_cov[k, l]
                    binwidth_i = hist_.GetBinWidth(i)
                    binwidth_j = hist_.GetBinWidth(j)
                    if i == k:
                        derivation_i = (integral - hist_.GetBinContent(i)) / pow(integral, 2) * (1 / binwidth_i)
                    else:
                        derivation_i = - hist_.GetBinContent(i) / pow(integral, 2) * (1 / binwidth_i)
                    if j == l:
                        derivation_j = (integral - hist_.GetBinContent(j)) / pow(integral, 2) * (1 / binwidth_j)
                    else:
                        derivation_j = - hist_.GetBinContent(j) / pow(integral, 2) * (1 / binwidth_j)
                    sum_value += derivation_i * derivation_j * old_entry
            new_cov[i, j] = sum_value

    return new_cov


def store_matrix_in_root(matrices, sample_names, pt_jet_ranges, filename):
    # type: (list, list, list, str) -> None

    f = ROOT.TFile(filename, 'RECREATE')
    f.cd()

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            for k, pt_jet_range in enumerate(pt_jet_ranges):
                hist = ROOT.TH2F("Covariance Matrix", "Normalized Covariance Matrix",
                                 matrices[g][h][k].shape[0], 0, 1, matrices[g][h][k].shape[1], 0, 1)
                for i in range(matrices[g][h][k].shape[0]):
                    for j in range(matrices[g][h][k].shape[1]):
                        hist.SetBinContent(i, j, matrices[g][h][k][i, j])

                hist.Write('norm_cov_matrix_{}_{}_{}-{}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1]))
                del hist

    f.Close()


pt_jet_lowest = 400
pt_jet_highest = 700
pt_jet_step = 50
pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                    range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))


if __name__ == '__main__':
    subfolder = '/cov_matrices'
    filename = 'histogram_files/correlator_hist_trip.root'
    sample_names = ['TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5']

    matrices = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            for k, pt_range in enumerate(pt_jet_ranges):
                matrices[g][h][k] = calc_norm_cov_matrix(filename_root_hist=filename,
                                                         hist_name='correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name,
                                                                                                            pt_range[0], pt_range[1]))

    # norm_cov_matrix = normalize_cov_matrix(matrix_orig=np.full((40, 40), 5), )
    store_matrix_in_root(matrices=matrices, sample_names=sample_names,
                         pt_jet_ranges=pt_jet_ranges, filename='cov_matrices/norm_cov_matrices.root')

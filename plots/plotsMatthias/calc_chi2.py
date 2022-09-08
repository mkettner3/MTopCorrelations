# -*- coding: utf-8 -*-

"""
Script to calculate the chi2 between histograms
"""

from typing import Any
import numpy as np
import ROOT
from MTopCorrelations.Tools.user import plot_directory


def calc_norm_cov_matrix(filename_root_hist, hist_name, plot_matrix=False, id_level=None, id_sample=None, id_range=None,
                         absolute_hist=False):
    # type: (str, str, bool, str, str, tuple, bool) -> tuple

    f = ROOT.TFile(filename_root_hist, 'read')
    root_hist = f.Get(hist_name)
    root_hist.SetDirectory(ROOT.nullptr)            # Returns a pointer to root_hist in memory.
    f.Close()                                       # f.Get only returns a handle, which gets lost when TFile is closed

    num_bins = root_hist.GetNbinsX()
    hist_axis_range_min = root_hist.GetXaxis().GetXmin()
    hist_axis_range_max = root_hist.GetXaxis().GetXmax()
    matrix_orig = np.zeros((num_bins, num_bins), dtype=np.float64)
    if absolute_hist:
        for i in range(num_bins):
            matrix_orig[i, i] = np.sqrt(root_hist.GetBinContent(i+1))
    else:
        for i in range(num_bins):
            matrix_orig[i, i] = root_hist.GetBinError(i+1) ** 2             # weight is set to Kronecker-Delta

    matrix_norm = normalize_cov_matrix(matrix_orig=matrix_orig, root_hist=root_hist)

    if plot_matrix:
        plot_matrix_in_root(matrix_norm, id_level, id_sample, id_range, (hist_axis_range_min, hist_axis_range_max),
                            absolute_hist)

    return matrix_norm, root_hist, (hist_axis_range_min, hist_axis_range_max)


def normalize_cov_matrix(matrix_orig, root_hist):
    # type: (np.ndarray, Any) -> np.ndarray

    num_bins = matrix_orig.shape[0]
    integral = root_hist.Integral()

    matrix_norm = np.zeros((num_bins, num_bins), dtype=np.float64)
    for i in range(num_bins):
        for j in range(num_bins):
            bincontent_i = root_hist.GetBinContent(i+1)
            bincontent_j = root_hist.GetBinContent(j+1)
            binwidth_i = root_hist.GetBinWidth(i+1)
            binwidth_j = root_hist.GetBinWidth(j+1)
            derivation_i_a = - bincontent_i / (integral**2 * binwidth_i)
            derivation_i_b = (integral - bincontent_i) / (integral**2 * binwidth_i)
            derivation_j_a = - bincontent_j / (integral**2 * binwidth_j)
            derivation_j_b = (integral - bincontent_j) / (integral**2 * binwidth_j)
            deriv_mat_i = np.full((num_bins, num_bins), derivation_i_a, dtype=np.float64)
            deriv_mat_i[i, :] = derivation_i_b
            deriv_mat_j = np.full((num_bins, num_bins), derivation_j_a, dtype=np.float64)
            deriv_mat_j[:, j] = derivation_j_b
            matrix_deriv = deriv_mat_i * deriv_mat_j * matrix_orig
            matrix_norm[i, j] = matrix_deriv.sum()

    return matrix_norm


def normalize_cov_matrix_python(matrix_orig, root_hist):
    # type: (np.ndarray, Any) -> np.ndarray

    num_bins = matrix_orig.shape[0]
    integral = root_hist.Integral()

    matrix_norm = np.zeros((num_bins, num_bins), dtype=np.float64)
    for i in range(num_bins):
        for j in range(num_bins):
            bincontent_i = root_hist.GetBinContent(i+1)
            bincontent_j = root_hist.GetBinContent(j+1)
            binwidth_i = root_hist.GetBinWidth(i+1)
            binwidth_j = root_hist.GetBinWidth(j+1)
            derivation_i_a = - bincontent_i / (integral**2 * binwidth_i)
            derivation_i_b = (integral - bincontent_i) / (integral**2 * binwidth_i)
            derivation_j_a = - bincontent_j / (integral**2 * binwidth_j)
            derivation_j_b = (integral - bincontent_j) / (integral**2 * binwidth_j)
            sum_element = 0
            for h in range(num_bins):
                for k in range(num_bins):
                    derivation_i = derivation_i_b if i == h else derivation_i_a
                    derivation_j = derivation_j_b if j == k else derivation_j_a
                    sum_element += derivation_i * derivation_j * matrix_orig[h, k]
            matrix_norm[i, j] = sum_element

    return matrix_norm


def NormalizeMatrix(old_cov, hist_):
    # type: (np.ndarray, Any) -> np.ndarray

    nbins = old_cov.shape[0]
    new_cov = np.empty((nbins, nbins), dtype=np.float64)
    integral = hist_.Integral()
    for i in range(nbins):
        for j in range(nbins):
            sum_value = 0
            for k in range(nbins):
                for l in range(nbins):
                    old_entry = old_cov[k, l]
                    binwidth_i = hist_.GetBinWidth(i+1)
                    binwidth_j = hist_.GetBinWidth(j+1)
                    if i == k:
                        derivation_i = (integral - hist_.GetBinContent(i+1)) / pow(integral, 2) * (1 / binwidth_i)
                    else:
                        derivation_i = - hist_.GetBinContent(i+1) / pow(integral, 2) * (1 / binwidth_i)
                    if j == l:
                        derivation_j = (integral - hist_.GetBinContent(j+1)) / pow(integral, 2) * (1 / binwidth_j)
                    else:
                        derivation_j = - hist_.GetBinContent(j+1) / pow(integral, 2) * (1 / binwidth_j)
                    sum_value += derivation_i * derivation_j * old_entry
            new_cov[i, j] = sum_value

    return new_cov


def compute_chi2(template_hist, data_hist, data_cov_matrix):
    # type: (Any, Any, np.ndarray) -> float

    num_bins = template_hist.GetNbinsX()
    template_hist_cont = np.zeros(num_bins, dtype=np.float64)
    for i in range(num_bins):
        template_hist_cont[i] = template_hist.GetBinContent(i+1)
    data_hist_cont = np.zeros(num_bins, dtype=np.float64)
    for i in range(num_bins):
        data_hist_cont[i] = data_hist.GetBinContent(i+1)
    d_vec = data_hist_cont - template_hist_cont

    data_cov_matrix_inv = np.linalg.inv(data_cov_matrix)
    chi2 = np.linalg.multi_dot([d_vec, data_cov_matrix_inv, d_vec])

    return chi2


def store_matrix_in_root(matrices_norm, matrices_orig, sample_names, pt_jet_ranges, filename, hist_axis_range):
    # type: (list, list, list, list, str, tuple) -> None

    f = ROOT.TFile(filename, 'RECREATE')
    f.cd()

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            for k, pt_jet_range in enumerate(pt_jet_ranges):
                hist_norm = ROOT.TH2F("Covariance Matrix", "Normalized Covariance Matrix",
                                      matrices_norm[g][h][k].shape[0], hist_axis_range[0], hist_axis_range[1],
                                      matrices_norm[g][h][k].shape[1], hist_axis_range[0], hist_axis_range[1])
                hist_orig = ROOT.TH2F("Covariance Matrix", "Original Covariance Matrix",
                                      matrices_orig[g][h][k].shape[0], hist_axis_range[0], hist_axis_range[1],
                                      matrices_orig[g][h][k].shape[1], hist_axis_range[0], hist_axis_range[1])
                for i in range(matrices_norm[g][h][k].shape[0]):
                    for j in range(matrices_norm[g][h][k].shape[1]):
                        hist_norm.SetBinContent(i, j, matrices_norm[g][h][k][i, j])
                for i in range(matrices_orig[g][h][k].shape[0]):
                    for j in range(matrices_orig[g][h][k].shape[1]):
                        hist_orig.SetBinContent(i, j, matrices_orig[g][h][k][i, j])

                hist_norm.Write('norm_cov_matrix_{}_{}_{}-{}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1]))
                hist_orig.Write('norm_cov_matrix_{}_{}_{}-{}_orig'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1]))
                del hist_norm
                del hist_orig

    f.Close()


def plot_matrix_in_root(matrix_norm, id_level, id_sample, id_range, hist_axis_range, absolute_hist):
    # type: (np.ndarray, str, str, tuple, tuple, bool) -> None

    hist_norm = ROOT.TH2F("Covariance Matrix", "Normalized Covariance Matrix",
                          matrix_norm.shape[0], hist_axis_range[0], hist_axis_range[1],
                          matrix_norm.shape[1], hist_axis_range[0], hist_axis_range[1])
    for i in range(matrix_norm.shape[0]):
        for j in range(matrix_norm.shape[1]):
            hist_norm.SetBinContent(i, j, matrix_norm[i, j])

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    c = ROOT.TCanvas('c', 'c', 1000, 1000)
    c.SetLogz()
    ROOT.gPad.SetRightMargin(0.2)
    hist_norm.SetTitle('Normalized Covariance Matrix')
    hist_norm.Draw('COLZ')

    abs = 'abs_' if absolute_hist else ''
    c.Print(plot_directory+'/corr_matrix_plots/correlation_matrix_norm_{}{:}_{:}_{:}-{:}.png'.format(abs, id_level, id_sample,
                                                                                                     id_range[0], id_range[1]))


def plot_chi2(root_graph, id_level, id_range):
    # type: (Any, str, tuple) -> None

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    c = ROOT.TCanvas('c', 'c', 1000, 1000)
    root_graph.SetTitle("Chi^{2}")
    root_graph.GetXaxis().SetTitle("Top-Mass (GeV)")
    # root_graph.GetYaxis().SetRangeUser(root_graph.GetMinimum()-1e10, root_graph.GetMaximum()+1e10)
    root_graph.Draw('')

    c.Print(plot_directory+'/chi2_plots/chi2_data_{}_{}-{}.png'.format(id_level, id_range[0], id_range[1]))


pt_jet_lowest = 400
pt_jet_highest = 700
pt_jet_step = 50
pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                    range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))


if __name__ == '__main__':
    subfolder = '/cov_matrices'
    filename = 'histogram_files/correlator_hist_trip_8.root'
    sample_names = ['TTbar_169p5', 'TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5', 'TTbar_175p5']

    matrices_norm = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    root_hist = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    chi2 = [[[] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    hist_axis_range = (0, 3)

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            print('Working on sample "{}" ...'.format(sample_name))
            for k, pt_range in enumerate(pt_jet_ranges):
                (matrices_norm[g][h][k],
                 root_hist[g][h][k],
                 hist_axis_range) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                         hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name,
                                                                                                            pt_range[0], pt_range[1]),
                                                         plot_matrix=False,
                                                         id_level=level, id_sample=sample_name, id_range=pt_range)

        for k, pt_range in enumerate(pt_jet_ranges):
            for h in [0, 1, 3, 4]:
                chi2[g][k].append(compute_chi2(template_hist=root_hist[g][h][k], data_hist=root_hist[g][2][k],
                                               data_cov_matrix=matrices_norm[g][2][k]))

            chi2_graph = ROOT.TGraph(4, np.array([169.5, 171.5, 173.5, 175.5]), np.asarray(chi2[g][k]))
            fit_func = ROOT.TF1('pol2_fit', 'pol2', 169.5, 175.5)
            fit_result = chi2_graph.Fit(fit_func, 'R')
            fit = chi2_graph.GetFunction('pol2_fit')
            plot_chi2(root_graph=chi2_graph, id_level=level, id_range=pt_range)
            obt_top_mass = fit.GetMinimumX()
            print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_mass))
            chi2min = fit.GetMinimum()
            uncertainty = fit.GetX(chi2min+1, 169.5, 175.5)
            print(uncertainty)

    """
    # f = ROOT.TFile(filename)
    # root_hist = f.Get('correlator_hist_Gen_TTbar_172p5_450_500')
    root_hist = ROOT.TH1F("Correlator", "3 #zeta", 5, 0, 3)
    root_hist.SetBinContent(1, 80)
    root_hist.SetBinContent(2, 40)
    root_hist.SetBinContent(3, 100)
    root_hist.SetBinContent(4, 10)
    root_hist.SetBinContent(5, 20)
    root_hist.SetBinError(1, 8)
    root_hist.SetBinError(2, 4)
    root_hist.SetBinError(3, 10)
    root_hist.SetBinError(4, 1)
    root_hist.SetBinError(5, 2)

    num_bins = root_hist.GetNbinsX()
    matrix_orig = np.zeros((num_bins, num_bins), dtype=np.float64)
    for i in range(num_bins):
        for j in range(num_bins):
            matrix_orig[i, j] = root_hist.GetBinError(i+1) * root_hist.GetBinError(j+1)

    norm_cov_matrix = normalize_cov_matrix(matrix_orig=matrix_orig, root_hist=root_hist)

    norm_cov_matrix = calc_norm_cov_matrix(filename_root_hist=filename, hist_name='correlator_hist_Gen_TTbar_172p5_450_500')
    print(norm_cov_matrix)
    """

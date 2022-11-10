# -*- coding: utf-8 -*-

"""
Script to calculate the chi2 between histograms with the option to scale the error value
"""

from typing import Any
from array import array
from copy import deepcopy
import numpy as np
import ROOT
from MTopCorrelations.Tools.user import plot_directory


def create_test_histograms(sample_name):
    if sample_name == 'TTbar_169p5':
        bin_con = [20, 80, 40, 10, 5]
    elif sample_name == 'TTbar_171p5':
        bin_con = [20, 70, 70, 10, 5]
    elif sample_name == 'TTbar_172p5':
        bin_con = [20, 50, 100, 10, 5]
    elif sample_name == 'TTbar_173p5':
        bin_con = [20, 40, 90, 50, 5]
    elif sample_name == 'TTbar_175p5':
        bin_con = [20, 30, 60, 80, 5]
    else:
        raise ValueError('The variable "sample_name" has an unexpected value!')

    root_hist = ROOT.TH1D("Correlator", "3 #zeta", 5, 0, 3)
    root_hist.SetBinContent(1, bin_con[0])
    root_hist.SetBinContent(2, bin_con[1])
    root_hist.SetBinContent(3, bin_con[2])
    root_hist.SetBinContent(4, bin_con[3])
    root_hist.SetBinContent(5, bin_con[4])
    root_hist.SetBinError(1, np.sqrt(bin_con[0]))
    root_hist.SetBinError(2, np.sqrt(bin_con[1]))
    root_hist.SetBinError(3, np.sqrt(bin_con[2]))
    root_hist.SetBinError(4, np.sqrt(bin_con[3]))
    root_hist.SetBinError(5, np.sqrt(bin_con[4]))

    return root_hist


def calc_norm_cov_matrix(filename_root_hist, hist_name, plot_matrix=False, id_level=None, id_sample=None, id_range=None,
                         absolute_hist=False, bin_error_scale=1):
    # type: (str, str, bool, str, str, tuple, bool, float) -> tuple

    histogram = prepare_histogram(filename_root_hist=filename_root_hist, hist_name=hist_name)

    num_bins = histogram.GetNbinsX()
    hist_axis_range_min = histogram.GetXaxis().GetXmin()
    hist_axis_range_max = histogram.GetXaxis().GetXmax()
    matrix_orig = np.zeros((num_bins, num_bins), dtype=np.float64)
    if absolute_hist:
        for i in range(num_bins):
            matrix_orig[i, i] = np.sqrt(histogram.GetBinContent(i+1))
    else:
        for i in range(num_bins):
            matrix_orig[i, i] = (histogram.GetBinError(i+1)*bin_error_scale) ** 2   # weight is set to Kronecker-Delta

    matrix_norm = normalize_cov_matrix(matrix_orig=matrix_orig, root_hist=histogram)

    if plot_matrix:
        plot_matrix_in_root(matrix_norm, id_level, id_sample, id_range, (hist_axis_range_min, hist_axis_range_max),
                            absolute_hist)

    return matrix_norm, histogram, (hist_axis_range_min, hist_axis_range_max)


def prepare_histogram(filename_root_hist, hist_name):
    # type: (str, str) -> Any

    f = ROOT.TFile(filename_root_hist, 'read')
    root_hist = f.Get(hist_name)
    root_hist.SetDirectory(ROOT.nullptr)            # Returns a pointer to root_hist in memory.
    f.Close()                                       # f.Get only returns a handle, which gets lost when TFile is closed

    hist_new = root_hist.Rebin(2, 'hist_new', array('d', [0.9, 2.03, 2.7]))

    return hist_new


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

    template_hist.Scale(1/template_hist.Integral(), 'width')
    data_hist.Scale(1/data_hist.Integral(), 'width')
    
    num_bins = template_hist.GetNbinsX()
    template_hist_cont = np.zeros((num_bins-1), dtype=np.float64)
    data_hist_cont = np.zeros((num_bins-1), dtype=np.float64)
    bin_list = list(range(num_bins))
    bin_list.remove(2)
    for idx, bin_nbr in enumerate(bin_list):
        template_hist_cont[idx] = template_hist.GetBinContent(bin_nbr+1)
        data_hist_cont[idx] = data_hist.GetBinContent(bin_nbr+1)
    d_vec = data_hist_cont - template_hist_cont

    data_cov_matrix = np.delete(np.delete(data_cov_matrix, 2, 0), 2, 1)
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
                hist_norm = ROOT.TH2D("Covariance Matrix", "Normalized Covariance Matrix",
                                      matrices_norm[g][h][k].shape[0], hist_axis_range[0], hist_axis_range[1],
                                      matrices_norm[g][h][k].shape[1], hist_axis_range[0], hist_axis_range[1])
                hist_orig = ROOT.TH2D("Covariance Matrix", "Original Covariance Matrix",
                                      matrices_orig[g][h][k].shape[0], hist_axis_range[0], hist_axis_range[1],
                                      matrices_orig[g][h][k].shape[1], hist_axis_range[0], hist_axis_range[1])
                for i in range(matrices_norm[g][h][k].shape[0]):
                    for j in range(matrices_norm[g][h][k].shape[1]):
                        hist_norm.SetBinContent(i+1, j+1, matrices_norm[g][h][k][i, j])
                for i in range(matrices_orig[g][h][k].shape[0]):
                    for j in range(matrices_orig[g][h][k].shape[1]):
                        hist_orig.SetBinContent(i+1, j+1, matrices_orig[g][h][k][i, j])

                hist_norm.Write('norm_cov_matrix_{}_{}_{}-{}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1]))
                hist_orig.Write('orig_cov_matrix_{}_{}_{}-{}'.format(level, sample_name, pt_jet_range[0], pt_jet_range[1]))
                del hist_norm
                del hist_orig

    f.Close()


def plot_matrix_in_root(matrix_norm, id_level, id_sample, id_range, hist_axis_range, absolute_hist):
    # type: (np.ndarray, str, str, tuple, tuple, bool) -> None

    hist_norm = ROOT.TH2D("Covariance Matrix", "Normalized Covariance Matrix",
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
    c.Print(plot_directory+'corr_matrix_plots/correlation_matrix_8_norm_{}{:}_{:}_{:}-{:}.png'.format(abs, id_level, id_sample,
                                                                                                      id_range[0], id_range[1]))


def plot_chi2(root_graph, label, filename, obt_top_mass):
    # type: (list, list, str, float) -> None

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    ROOT.gStyle.SetLegendBorderSize(1)  # No border for legend
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(0)
    c = ROOT.TCanvas('c', 'c', 1000, 1000)
    legend = ROOT.TLegend(0.3, 0.74, 0.7, 0.9)

    graphs = ROOT.TMultiGraph()
    for s in range(len(root_graph)):
        root_graph[s].SetMarkerSize(2)
        root_graph[s].SetMarkerStyle(47)
        root_graph[s].SetMarkerColor(s+2)
        root_graph[s].GetFunction('pol2_fit').SetLineColor(s+2)
        graphs.Add(deepcopy(root_graph[s]))
        legend.AddEntry(root_graph[s].GetFunction('pol2_fit'), label[s], 'l')
    graphs.SetTitle('#chi^{2}')
    graphs.GetXaxis().SetTitle('Top-Mass (GeV)')
    graphs.SetMaximum(max([root_graph[i].GetHistogram().GetMaximum() for i in range(len(root_graph))])+20)
    graphs.SetMinimum(min([root_graph[j].GetHistogram().GetMinimum() for j in range(len(root_graph))])-20)
    graphs.Draw('AP')
    legend.AddEntry(ROOT.nullptr, 'Resulting Top-Mass: {:.3f} GeV'.format(obt_top_mass), '')
    legend.Draw()

    c.Print(plot_directory+filename)


def plot_corr_hist(corr_hists, hist_range, filename_graphic, sample_names):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    c = ROOT.TCanvas('c', 'c', 600, 600)
    legend = ROOT.TLegend(0.75, 0.7, 0.94, 0.89)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    if len(corr_hists) == 3:
        line_colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kRed]
    elif len(corr_hists) == 5:
        line_colors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kYellow]
    else:
        raise RuntimeError('Please specify the line colors in plot_corr_hist()!')

    for hist, line_color, sample_name in zip(corr_hists, line_colors, sample_names):
        hist.SetLineColor(line_color)
        hist.SetTitle('')
        hist.SetLineWidth(2)
        hist.SetLineStyle(1)
        legend.AddEntry(hist, sample_name[-5:], 'l')
    corr_hists[0].GetXaxis().SetRangeUser(hist_range[0], hist_range[1])
    corr_hists[0].GetXaxis().SetTitle('3#zeta')
    corr_hists[0].GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    corr_hists[0].GetYaxis().SetRangeUser(0, max([corr_hists[i].GetMaximum() for i in range(len(corr_hists))])*1.1)
    corr_hists[0].GetYaxis().SetTitle('Energy-weighted Triplets')
    corr_hists[0].GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse

    corr_hists[0].Draw('HIST')
    for i in range(1, len(corr_hists)):
        corr_hists[i].Draw('HIST SAME')

    legend.Draw()
    c.Print(plot_directory+filename_graphic)


pt_jet_lowest = 400
pt_jet_highest = 700
pt_jet_step = 50
pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                    range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))


if __name__ == '__main__':
    filename = 'histogram_files/correlator_hist_trip_16.root'
    sample_names = ['TTbar_169p5', 'TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5', 'TTbar_175p5']
    error_scales = [1, 2, 4]

    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    matrices_norm = [[[[None for _ in range(len(error_scales))] for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    root_hist = [[[[None for _ in range(len(error_scales))] for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    chi2 = [[[[] for _ in range(len(error_scales))] for _ in range(len(pt_jet_ranges))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            print('Working on sample "{}" ...'.format(sample_name))
            for k, pt_range in enumerate(pt_jet_ranges):
                for s, scale in enumerate(error_scales):
                    (matrices_norm[g][h][k][s],
                     root_hist[g][h][k][s],
                     hist_range) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                        hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name,
                                                                                        pt_range[0], pt_range[1]),
                                                        plot_matrix=False,
                                                        id_level=level, id_sample=sample_name, id_range=pt_range,
                                                        bin_error_scale=scale)

        for k, pt_range in enumerate(pt_jet_ranges):
            plot_corr_hist(corr_hists=[root_hist[g][i][k][0] for i in range(len(sample_names))], hist_range=hist_range,
                           filename_graphic='chi2_plots/chi2_new_16_hist/corr_hist_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                           sample_names=sample_names)
            for s in range(len(error_scales)):
                for h in [0, 1, 3, 4]:
                    chi2[g][k][s].append(compute_chi2(template_hist=root_hist[g][h][k][s], data_hist=root_hist[g][2][k][s],
                                                      data_cov_matrix=matrices_norm[g][2][k][s]))

            chi2_graph = [ROOT.TGraph(4, np.array([169.5, 171.5, 173.5, 175.5]), np.asarray(chi2[g][k][s])) for s in range(len(error_scales))]
            fit_func = ROOT.TF1('pol2_fit', 'pol2', 168, 177)
            [chi2_graph[s].Fit(fit_func, 'R') for s in range(len(error_scales))]
            fit = [chi2_graph[s].GetFunction('pol2_fit') for s in range(len(error_scales))]
            obt_top_mass = fit[0].GetMinimumX()
            print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_mass))
            chi2min = fit[0].GetMinimum()
            uncertainty = abs(obt_top_mass - fit[0].GetX(chi2min+1, 168, 177))
            print('The uncertainty equals {:.5f} GeV.'.format(uncertainty))
            plot_chi2(root_graph=chi2_graph, label=['Error factor: {:0}'.format(e) for e in error_scales],
                      filename='chi2_plots/chi2_new_16_{}_{}-{}.pdf'.format(level, pt_range[0], pt_range[1]),
                      obt_top_mass=obt_top_mass)

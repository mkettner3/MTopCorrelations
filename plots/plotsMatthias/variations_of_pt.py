# -*- coding: utf-8 -*-

"""
Script to determine the systematic error caused by variations of pT of a jet and its constituents.
"""

from calc_triplet_and_hist import pt_jet_ranges
from calc_chi2 import calc_norm_cov_matrix, normalize_cov_matrix, compute_chi2
import numpy as np
import ROOT

# subfolder = '/cov_matrices'
filename = 'histogram_files/correlator_hist_trip_9.root'
sample_names = ['TTbar_169p5', 'TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5', 'TTbar_175p5']

matrices_norm = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
root_hist = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
hists_varied = [[[None for _ in range(2)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
sigma_up = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
sigma_down = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
matrices_norm_up = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
matrices_norm_down = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
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
    f = ROOT.TFile(filename, 'read')
    for k, pt_jet_range in enumerate(pt_jet_ranges):
        for v, var_fac in enumerate(['1.02', '0.98']):
            hists_varied[g][k][v] = f.Get('/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_{:}_{:}_{:}_{:}'.format(var_fac, level, pt_jet_range[0], pt_jet_range[1]))
            hists_varied[g][k][v].SetDirectory(ROOT.nullptr)
    f.Close()

    num_bins = root_hist[g][0][0].GetNbinsX()
    hist_axis_range_min = root_hist[g][0][0].GetXaxis().GetXmin()
    hist_axis_range_max = root_hist[g][0][0].GetXaxis().GetXmax()
    for k, pt_jet_range in enumerate(pt_jet_ranges):
        sigma_up[g][k] = [hists_varied[g][k][0].GetBinContent(i+1) - root_hist[g][2][k].GetBinContent(i+1) for i in range(num_bins)]
        sigma_down[g][k] = [hists_varied[g][k][1].GetBinContent(i+1) - root_hist[g][2][k].GetBinContent(i+1) for i in range(num_bins)]

        matrix_orig_up = np.zeros((num_bins, num_bins), dtype=np.float64)
        matrix_orig_down = np.zeros((num_bins, num_bins), dtype=np.float64)
        for i in range(num_bins):
            for j in range(num_bins):
                matrix_orig_up[i, j] = sigma_up[g][k][i] * sigma_up[g][k][j]
                matrix_orig_down[i, j] = sigma_down[g][k][i] * sigma_down[g][k][j]

        matrices_norm_up[g][k] = normalize_cov_matrix(matrix_orig=matrix_orig_up, root_hist=hists_varied[g][k][0]) + matrices_norm[g][2][k]
        matrices_norm_down[g][k] = normalize_cov_matrix(matrix_orig=matrix_orig_down, root_hist=hists_varied[g][k][1])  + matrices_norm[g][2][k]

    uncertainty_stat = []
    uncertainty_tot = []
    for k, pt_range in enumerate(pt_jet_ranges):
        for h in [0, 1, 3, 4]:
            chi2[g][k].append(compute_chi2(template_hist=root_hist[g][h][k], data_hist=root_hist[g][2][k],
                                           data_cov_matrix=matrices_norm[g][2][k]))

        chi2_graph = ROOT.TGraph(4, np.array([169.5, 171.5, 173.5, 175.5]), np.asarray(chi2[g][k]))
        fit_func = ROOT.TF1('pol2_fit', 'pol2', 169.5, 175.5)
        fit_result = chi2_graph.Fit(fit_func, 'R')
        fit = chi2_graph.GetFunction('pol2_fit')
        # plot_chi2(root_graph=chi2_graph, id_level=level, id_range=pt_range)
        obt_top_mass = fit.GetMinimumX()
        print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_mass))
        chi2min = fit.GetMinimum()
        uncertainty_stat.append(fit.GetX(chi2min+1, 169.5, 175.5))
        print(uncertainty_stat[-1])

    for k, pt_range in enumerate(pt_jet_ranges):
        for h in [0, 1, 3, 4]:
            chi2[g][k].append(compute_chi2(template_hist=root_hist[g][h][k], data_hist=hists_varied[g][k][0],
                                           data_cov_matrix=matrices_norm_up[g][k]))

        chi2_graph = ROOT.TGraph(4, np.array([169.5, 171.5, 173.5, 175.5]), np.asarray(chi2[g][k]))
        fit_func = ROOT.TF1('pol2_fit', 'pol2', 169.5, 175.5)
        fit_result = chi2_graph.Fit(fit_func, 'R')
        fit = chi2_graph.GetFunction('pol2_fit')
        # plot_chi2(root_graph=chi2_graph, id_level=level, id_range=pt_range)
        obt_top_mass = fit.GetMinimumX()
        print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_mass))
        chi2min = fit.GetMinimum()
        uncertainty_tot.append(fit.GetX(chi2min+1, 169.5, 175.5))
        print(uncertainty_tot[-1])

    uncertainty_up = [np.sqrt(uncertainty_tot[k]**2 - uncertainty_stat[k]**2) for k in range(len(pt_jet_ranges))]
    print('Uncertainty Up:')
    for k in range(len(pt_jet_ranges)):
        print(uncertainty_up[k])

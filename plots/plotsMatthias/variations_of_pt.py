# -*- coding: utf-8 -*-

"""
Script to determine the systematic error caused by variations of pT of a jet and its constituents.
"""

from calc_chi2 import pt_jet_ranges, calc_norm_cov_matrix, prepare_histogram, normalize_cov_matrix, compute_chi2,\
    plot_chi2, plot_corr_hist
import numpy as np
import ROOT

if __name__ == '__main__':
    filename = 'histogram_files/correlator_hist_trip_13.root'
    sample_names = ['TTbar_169p5', 'TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5', 'TTbar_175p5']

    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    matrices_norm = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    root_hist = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    hists_varied = [[[None for _ in range(2)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    sigma_up = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    sigma_down = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    matrices_norm_up = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    matrices_norm_down = [[None for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    chi2 = [[[[] for _ in range(3)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            print('Working on sample "{}" ...'.format(sample_name))
            for k, pt_range in enumerate(pt_jet_ranges):
                (matrices_norm[g][h][k],
                 root_hist[g][h][k],
                 hist_range) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                    hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name, pt_range[0], pt_range[1]),
                                                    plot_matrix=False,
                                                    id_level=level, id_sample=sample_name, id_range=pt_range)

        for k, pt_jet_range in enumerate(pt_jet_ranges):
            for v, var_fac in enumerate((1.02, 0.98)):
                hists_varied[g][k][v] = prepare_histogram(filename_root_hist=filename,
                                                          hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_{:.2f}_{:}_{:}_{:}_{:}'.format(var_fac, level, 'TTbar_172p5', pt_jet_range[0], pt_jet_range[1]))

        num_bins = root_hist[g][0][0].GetNbinsX()
        uncertainty_stat = []
        uncertainty_tot = []
        for k, pt_jet_range in enumerate(pt_jet_ranges):
            plot_corr_hist(corr_hists=[root_hist[g][2][k], hists_varied[g][k][0], hists_varied[g][k][1]], hist_range=hist_range,
                           filename_graphic='chi2_plots/chi2_pt_varied_13_hist/corr_hist_{}_{}-{}.png'.format(level, pt_jet_range[0], pt_jet_range[1]),
                           sample_names=['p_{T} variance: '+e for e in ['original', '+ 2 %', '- 2 %']])
            sigma_up[g][k] = [hists_varied[g][k][0].GetBinContent(i+1) - root_hist[g][2][k].GetBinContent(i+1) for i in range(num_bins)]
            sigma_down[g][k] = [hists_varied[g][k][1].GetBinContent(i+1) - root_hist[g][2][k].GetBinContent(i+1) for i in range(num_bins)]

            matrix_orig_up = np.zeros((num_bins, num_bins), dtype=np.float64)
            matrix_orig_down = np.zeros((num_bins, num_bins), dtype=np.float64)
            for i in range(num_bins):
                for j in range(num_bins):
                    matrix_orig_up[i, j] = sigma_up[g][k][i] * sigma_up[g][k][j]
                    matrix_orig_down[i, j] = sigma_down[g][k][i] * sigma_down[g][k][j]

            matrices_norm_up[g][k] = normalize_cov_matrix(matrix_orig=matrix_orig_up, root_hist=root_hist[g][2][k]) + matrices_norm[g][2][k]
            matrices_norm_down[g][k] = normalize_cov_matrix(matrix_orig=matrix_orig_down, root_hist=root_hist[g][2][k]) + matrices_norm[g][2][k]

            for h in [0, 1, 3, 4]:
                chi2[g][k][0].append(compute_chi2(template_hist=root_hist[g][h][k], data_hist=root_hist[g][2][k],
                                                  data_cov_matrix=matrices_norm[g][2][k]))
                chi2[g][k][1].append(compute_chi2(template_hist=root_hist[g][h][k], data_hist=root_hist[g][2][k],
                                                  data_cov_matrix=matrices_norm_up[g][k]))
                chi2[g][k][2].append(compute_chi2(template_hist=root_hist[g][h][k], data_hist=root_hist[g][2][k],
                                                  data_cov_matrix=matrices_norm_down[g][k]))

            chi2_graph = [ROOT.TGraph(4, np.array([169.5, 171.5, 173.5, 175.5]), np.asarray(chi2[g][k][v])) for v in range(3)]
            fit_func = ROOT.TF1('pol2_fit', 'pol2', 169.5, 175.5)
            [chi2_graph[v].Fit(fit_func, 'R') for v in range(3)]
            fit = [chi2_graph[v].GetFunction('pol2_fit') for v in range(3)]
            obt_top_masses = [fit[v].GetMinimumX() for v in range(3)]
            print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_masses[0]))
            chi2min = [fit[v].GetMinimum() for v in range(3)]
            uncertainty_stat = abs(obt_top_masses[0] - fit[0].GetX(chi2min[0]+1, 169.5, 175.5))
            uncertainty_tot = [abs(obt_top_masses[v] - fit[v].GetX(chi2min[v]+1, 169.5, 175.5)) for v in range(1, 3)]
            print('The total uncertainty equals {:.5f} GeV.'.format(uncertainty_tot[0]))

            uncertainty_up = np.sqrt(uncertainty_tot[0]**2 - uncertainty_stat**2)
            uncertainty_down = np.sqrt(uncertainty_tot[1]**2 - uncertainty_stat**2)
            print('Uncertainty Up: {:.5f} GeV.'.format(uncertainty_up))
            print('Uncertainty Down: {:.5f} GeV.'.format(uncertainty_down))

            plot_chi2(root_graph=chi2_graph, label=['p_{T} variance: '+e for e in ['original', '+ 2 %', '- 2 %']],
                      filename='chi2_plots/chi2_pt_varied_13_{}_{}-{}.pdf'.format(level, pt_jet_range[0], pt_jet_range[1]),
                      obt_top_mass=obt_top_masses[1])

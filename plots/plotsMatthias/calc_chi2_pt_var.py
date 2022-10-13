# -*- coding: utf-8 -*-

"""
Script to calculate the chi2 between histograms with the option to scale the error value
"""

from calc_chi2 import pt_jet_ranges, calc_norm_cov_matrix, compute_chi2, plot_chi2, plot_corr_hist
import numpy as np
import ROOT


if __name__ == '__main__':
    filename = 'histogram_files/correlator_hist_trip_10.root'
    sample_names = ['TTbar_169p5', 'TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5', 'TTbar_175p5']

    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    matrices_norm = [[[[None for _ in range(3)] for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    root_hist = [[[[None for _ in range(3)] for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    chi2 = [[[[] for _ in range(3)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            print('Working on sample "{}" ...'.format(sample_name))
            for k, pt_range in enumerate(pt_jet_ranges):
                (matrices_norm[g][h][k][0],
                 root_hist[g][h][k][0]) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                               hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name,
                                                                                                                pt_range[0], pt_range[1]),
                                                               plot_matrix=False,
                                                               id_level=level, id_sample=sample_name, id_range=pt_range)
                for v, var_fac in enumerate((1.02, 0.98)):
                    (matrices_norm[g][h][k][v+1],
                     root_hist[g][h][k][v+1]) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                                     hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_{:.2f}_{:}_{:}_{:}'.format(var_fac, level, pt_range[0], pt_range[1]),
                                                                     plot_matrix=False,
                                                                     id_level=level, id_sample=sample_name, id_range=pt_range)

        for k, pt_range in enumerate(pt_jet_ranges):
            plot_corr_hist(corr_hists=[root_hist[g][i][k][0] for i in range(len(sample_names))],
                           filename_graphic='chi2_plots/chi2_pt_varied_10_hist/corr_hist_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]))
            for v in range(3):
                for h in [0, 1, 3, 4]:
                    chi2[g][k][v].append(compute_chi2(template_hist=root_hist[g][h][k][0], data_hist=root_hist[g][2][k][v],         # Muss v für template_hist 0 gesetzt werden oder nicht?
                                                      data_cov_matrix=matrices_norm[g][2][k][v]))

            chi2_graph = [ROOT.TGraph(4, np.array([169.5, 171.5, 173.5, 175.5]), np.asarray(chi2[g][k][v])) for v in range(3)]
            fit_func = ROOT.TF1('pol2_fit', 'pol2', 169.5, 175.5)
            [chi2_graph[v].Fit(fit_func, 'R') for v in range(3)]
            fit = [chi2_graph[v].GetFunction('pol2_fit') for v in range(3)]
            obt_top_mass = fit[0].GetMinimumX()
            print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_mass))
            chi2min = fit[0].GetMinimum()
            uncertainty = fit[0].GetX(chi2min+1, 169.5, 175.5)
            print(uncertainty)
            plot_chi2(root_graph=chi2_graph, label=['p_{T} variance: '+e for e in ['original', '+ 2 %', '- 2 %']],
                      filename='chi2_plots/chi2_pt_varied_n_10_{}_{}-{}.pdf'.format(level, pt_range[0], pt_range[1]),
                      obt_top_mass=obt_top_mass)

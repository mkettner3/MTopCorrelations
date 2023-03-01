# -*- coding: utf-8 -*-

"""
Script to determine the systematic error caused by variations of pT of a jet and its constituents.
"""

from calc_chi2 import pt_jet_ranges, calc_norm_cov_matrix, prepare_histogram, normalize_cov_matrix, compute_chi2,\
    plot_corr_hist
from MTopCorrelations.Tools.user import plot_directory
from copy import deepcopy
import numpy as np
import ROOT


def plot_chi2(root_graph, label, filename, obt_top_masses, uncertainties):
    # type: (list, list, str, list, list) -> None

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    ROOT.gStyle.SetLegendBorderSize(1)  # No border for legend
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(0)
    c = ROOT.TCanvas('c', 'c', 1000, 1000)
    legend = ROOT.TLegend(0.2, 0.6, 0.8, 0.9)

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
    legend.AddEntry(ROOT.nullptr, 'Resulting Top-Mass (+2%): {:.3f} GeV'.format(obt_top_masses[0]), '')
    legend.AddEntry(ROOT.nullptr, 'Uncertainty (+2%): {:.3f} GeV'.format(uncertainties[0]), '')
    legend.AddEntry(ROOT.nullptr, 'Resulting Top-Mass (-2%): {:.3f} GeV'.format(obt_top_masses[1]), '')
    legend.AddEntry(ROOT.nullptr, 'Uncertainty (-2%): {:.3f} GeV'.format(uncertainties[1]), '')
    legend.Draw()

    c.Print(plot_directory+filename)


def plot_uncertainties(uncertainties, filename):
    # type: (list, str) -> None

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    ROOT.gStyle.SetLegendBorderSize(1)  # No border for legend
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(0)
    c = ROOT.TCanvas('c', 'c', 1000, 1000)

    uncertainties = [np.mean([uncertainties[i], uncertainties[-(i+1)]]) for i in range(int(len(uncertainties)/2))]

    graph = ROOT.TGraph(4, np.array([1.1, 1.05, 1.02, 1.01]), np.asarray(uncertainties))
    graph.SetMarkerSize(2)
    graph.SetMarkerStyle(47)
    graph.SetMarkerColor(2)
    graph.SetLineColor(2)
    graph.SetTitle('Uncertainties for varied p_{T}')
    graph.GetXaxis().SetTitle('Variation of p_{T}')
    graph.Draw('APL')

    c.Print(plot_directory+filename)


if __name__ == '__main__':
    filename = 'histogram_files/correlator_hist_trip_26.root'
    sample_names = ['171.5', '171.75', '172.0', '172.25', 'None', '172.75', '173.0', '173.25', '173.5']

    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    matrices_norm = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    root_hist = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    root_hist_norm = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
    hists_varied = [[[None for _ in range(8)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    hists_varied_norm = [[[None for _ in range(8)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    sigma = [[[None for _ in range(8)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    matrices_varied_norm = [[[None for _ in range(8)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]
    chi2 = [[[None for _ in range(9)] for _ in range(len(pt_jet_ranges))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample_name in enumerate(sample_names):
            print('Working on sample "{}" ...'.format(sample_name))
            for k, pt_range in enumerate(pt_jet_ranges):
                (matrices_norm[g][h][k],
                 root_hist[g][h][k],
                 root_hist_norm[g][h][k],
                 hist_range) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                    hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name, pt_range[0], pt_range[1]),
                                                    plot_matrix=False,
                                                    id_level=level, id_sample=sample_name, id_range=pt_range)

        for k, pt_jet_range in enumerate(pt_jet_ranges):
            for v, var_fac in enumerate([1.1, 1.05, 1.02, 1.01, 0.99, 0.98, 0.95, 0.9]):
                hists_varied[g][k][v] = prepare_histogram(filename_root_hist=filename,
                                                          hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_jet_{:.2f}_{:}_{:}_{:}_{:}'.format(var_fac, level, 'None', pt_jet_range[0], pt_jet_range[1]))
                hists_varied_norm[g][k][v] = hists_varied[g][k][v].Clone()
                hists_varied_norm[g][k][v].Scale(1 / hists_varied_norm[g][k][v].Integral(), 'width')

        num_bins = root_hist[g][0][0].GetNbinsX()
        for k, pt_jet_range in enumerate(pt_jet_ranges):
            for v in range(8):
                sigma[g][k][v] = [hists_varied[g][k][v].GetBinContent(i+1) - root_hist[g][4][k].GetBinContent(i+1) for i in range(num_bins)]

                matrix_varied_orig = [np.zeros((num_bins, num_bins), dtype=np.float64) for _ in range(8)]
                for i in range(num_bins):
                    for j in range(num_bins):
                        matrix_varied_orig[v][i, j] = sigma[g][k][v][i] * sigma[g][k][v][j]

                matrices_varied_norm[g][k][v] = normalize_cov_matrix(matrix_orig=matrix_varied_orig[v], root_hist=root_hist[g][4][k]) + matrices_norm[g][4][k]

                plot_corr_hist(corr_hists=[root_hist_norm[g][4][k], hists_varied_norm[g][k][2], hists_varied_norm[g][k][5]], hist_range=hist_range,
                               filename_graphic='chi2_plots/chi2_pt_varied_26_hist/corr_hist_{}_{}-{}.png'.format(level, pt_jet_range[0], pt_jet_range[1]),
                               sample_names=['p_{T} variance: '+e for e in ['original', '+ 2 %', '- 2 %']])

            chi2[g][k][0] = [compute_chi2(template_hist=root_hist_norm[g][h][k], data_hist=root_hist_norm[g][4][k], data_cov_matrix=matrices_norm[g][4][k]) for h in range(9)]
            for v in range(8):
                chi2[g][k][v+1] = [compute_chi2(template_hist=root_hist_norm[g][h][k], data_hist=root_hist_norm[g][4][k], data_cov_matrix=matrices_varied_norm[g][k][v]) for h in range(9)]

            chi2_graph = [ROOT.TGraph(9, np.array([171.5, 171.75, 172.0, 172.25, 172.5, 172.75, 173.0, 173.25, 173.5]), np.asarray(chi2[g][k][v])) for v in range(9)]
            fit_func = ROOT.TF1('pol2_fit', 'pol2', 170, 175)
            [chi2_graph[v].Fit(fit_func, 'R') for v in range(9)]
            fit = [chi2_graph[v].GetFunction('pol2_fit') for v in range(9)]
            obt_top_masses = [fit[v].GetMinimumX() for v in range(9)]
            print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_masses[0]))
            chi2min = [fit[v].GetMinimum() for v in range(9)]
            uncertainty_stat = abs(obt_top_masses[0] - fit[0].GetX(chi2min[0]+1, 170, 175))
            uncertainty_tot = [abs(obt_top_masses[v] - fit[v].GetX(chi2min[v]+1, 170, 175)) for v in range(1,9)]
            print('The total uncertainty equals {:.5f} GeV.'.format(uncertainty_tot[2]))

            uncertainties = [np.sqrt(uncertainty_tot[v]**2 - uncertainty_stat**2) for v in range(8)]
            print('Uncertainty Up: {:.5f} GeV.'.format(uncertainties[2]))
            print('Uncertainty Down: {:.5f} GeV.'.format(uncertainties[5]))

            plot_chi2(root_graph=[chi2_graph[i] for i in range(9)], label=['p_{T} variance: '+e for e in ['original', '+ 10 %', '+ 5 %', '+ 2 %', '+ 1 %', '- 1 %', '- 2 %', '- 5 %', '- 10 %']],
                      filename='chi2_plots/chi2_pt_varied_26/chi2_pt_varied_26_{}_{}-{}.pdf'.format(level, pt_jet_range[0], pt_jet_range[1]),
                      obt_top_masses=[obt_top_masses[3], obt_top_masses[6]], uncertainties=[uncertainties[2], uncertainties[5]])

            plot_uncertainties(uncertainties=uncertainties,
                               filename='chi2_plots/chi2_pt_varied_26/chi2_pt_varied_26_uncertainties_{}_{}-{}.pdf'.format(level, pt_jet_range[0], pt_jet_range[1]))

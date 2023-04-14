# -*- coding: utf-8 -*-

"""
Script to determine the systematic error caused by variations of pT of a jet and its constituents.
"""

from calc_chi2 import pt_jet_ranges, prepare_histogram, normalize_cov_matrix, compute_chi2,\
    plot_corr_hist, plot_matrix_in_root, plot_vector_in_root
from MTopCorrelations.Tools.user import plot_directory
from copy import deepcopy
import numpy as np
import ROOT


def plot_constituent_hist(filename_root, histogram_name, hist_range, filename_graphic):
    f = ROOT.TFile(filename_root, 'read')
    histogram = f.Get(histogram_name)
    histogram.SetDirectory(ROOT.nullptr)            # Returns a pointer to root_hist in memory.
    f.Close()                                       # f.Get only returns a handle, which gets lost when TFile is closed

    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    histogram.Rebin(20)
    histogram.GetXaxis().SetRangeUser(hist_range[0], hist_range[1])
    histogram.GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    histogram.GetXaxis().SetTitle('p_{T} (GeV)')
    histogram.GetYaxis().SetRangeUser(0, histogram.GetMaximum()*1.1)
    histogram.GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    histogram.SetTitle('Constituent p_{T}')

    histogram.Draw('HIST')
    c.Print(plot_directory+filename_graphic)


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
    if obt_top_masses is not None or uncertainties is not None:
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

    # uncertainties = [np.mean([uncertainties[i], uncertainties[-(i+1)]]) for i in range(int(len(uncertainties)/2))]

    graph = ROOT.TGraph(len(var_factors), np.array(var_factors), np.asarray(uncertainties))
    graph.SetMarkerSize(2)
    graph.SetMarkerStyle(47)
    graph.SetMarkerColor(2)
    graph.SetLineColor(2)
    graph.SetTitle('Uncertainties for varied p_{T}')
    graph.GetXaxis().SetTitle('Variation of p_{T}')
    graph.Draw('APL')

    c.Print(plot_directory+filename)


def main(level, pt_range):
    corr_hist = [prepare_histogram(filename_root_hist=filename,
                                   hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name, pt_range[0], pt_range[1]),
                                   hist_binning=hist_binning) for sample_name in sample_names_stored]

    corr_hist_varied_jet = [prepare_histogram(filename_root_hist=filename,
                                              hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_jet_{:.2f}_{:}_None_{:}_{:}'.format(var_fac, level,
                                                                                                                                                        pt_range[0],
                                                                                                                                                        pt_range[1]),
                                              hist_binning=hist_binning) for var_fac in var_factors]

    corr_hist_varied_cons = [prepare_histogram(filename_root_hist=filename,
                                               hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_cons_pt_{:.2f}_{:}_None_{:}_{:}'.format(var_fac, level,
                                                                                                                                                             pt_range[0],
                                                                                                                                                             pt_range[1]),
                                               hist_binning=hist_binning) for var_fac in var_factors]

    plot_corr_hist(corr_hists=corr_hist,
                   filename_graphic='chi2_plots/chi2_pt_varied_28_hist/corr_hist_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                   sample_names=sample_names, title='Correlator Histograms (BW-reweighted)', hist_range=None)
    plot_corr_hist(corr_hists=[corr_hist[4]]+[corr_hist_varied_jet[v] for v in range(len(var_factors))],
                   filename_graphic='chi2_plots/chi2_pt_varied_28_hist/corr_hist_varied_jet_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                   sample_names=['original']+var_factor_names, title='Jet-p_{T} Variances',  hist_range=None)
    plot_corr_hist(corr_hists=[corr_hist[4]]+[corr_hist_varied_cons[v] for v in range(len(var_factors))],
                   filename_graphic='chi2_plots/chi2_pt_varied_28_hist/corr_hist_varied_cons_pt_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                   sample_names=['original']+var_factor_names, title='Constituent-p_{T} Variances',  hist_range=None)

    num_bins = corr_hist[0].GetNbinsX()
    matrix_orig = [np.zeros((num_bins, num_bins), dtype=np.float64) for _ in range(len(sample_names))]
    for h in range(len(sample_names)):
        for i in range(num_bins):
            matrix_orig[h][i, i] = (corr_hist[h].GetBinError(i+1)) ** 2   # weight is set to Kronecker-Delta

    sigma_jet = [[corr_hist_varied_jet[v].GetBinContent(i+1) - corr_hist[4].GetBinContent(i+1) for i in range(num_bins)] for v in range(len(var_factors))]
    sigma_cons = [[corr_hist_varied_cons[v].GetBinContent(i+1) - corr_hist[4].GetBinContent(i+1) for i in range(num_bins)] for v in range(len(var_factors))]

    matrix_varied_jet_orig = [np.zeros((num_bins, num_bins), dtype=np.float64) for _ in range(len(var_factors))]
    for v in range(len(var_factors)):
        for i in range(num_bins):
            for j in range(num_bins):
                matrix_varied_jet_orig[v][i, j] = sigma_jet[v][i] * sigma_jet[v][j]

    matrix_varied_cons_orig = [np.zeros((num_bins, num_bins), dtype=np.float64) for _ in range(len(var_factors))]
    for v in range(len(var_factors)):
        for i in range(num_bins):
            for j in range(num_bins):
                matrix_varied_cons_orig[v][i, j] = sigma_cons[v][i] * sigma_cons[v][j]

    for h, sample_name in enumerate(sample_names):
        plot_matrix_in_root(matrix=matrix_orig[h], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/matrix_orig_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], sample_name),
                            hist_axis_range=hist_binning, title='Covariance Matrix Sample '+sample_name)

    for v, var_name in enumerate(var_factor_names):
        plot_vector_in_root(vector=sigma_jet[v], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/corr_sigmas_jet_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], var_name),
                            hist_axis_range=hist_binning, title='Correlator Sigmas Jet-p_{T}-variation '+var_name)

    for v, var_name in enumerate(var_factor_names):
        plot_matrix_in_root(matrix=matrix_varied_jet_orig[v], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/matrix_varied_jet_orig_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], var_name),
                            hist_axis_range=hist_binning, title='Covariance Matrix Jet-p_{T}-variation '+var_name)

    for v, var_name in enumerate(var_factor_names):
        plot_vector_in_root(vector=sigma_cons[v], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/corr_sigmas_cons_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], var_name),
                            hist_axis_range=hist_binning, title='Correlator Sigmas Constituent-p_{T}-variation '+var_name)

    for v, var_name in enumerate(var_factor_names):
        plot_matrix_in_root(matrix=matrix_varied_cons_orig[v], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/matrix_varied_cons_orig_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], var_name),
                            hist_axis_range=hist_binning, title='Covariance Matrix Constituent-p_{T}-variation '+var_name)

    corr_hist_norm = [corr_hist[h].Clone() for h in range(len(sample_names))]
    for h in range(len(sample_names)):
        corr_hist_norm[h].Scale(1 / corr_hist_norm[h].Integral(), 'width')

    matrix_norm = [normalize_cov_matrix(matrix_orig=matrix_orig[h], root_hist=corr_hist[4]) for h in range(len(sample_names))]

    corr_hist_varied_jet_norm = [corr_hist_varied_jet[v].Clone() for v in range(len(var_factors))]
    for v in range(len(var_factors)):
        corr_hist_varied_jet_norm[v].Scale(1 / corr_hist_varied_jet_norm[v].Integral(), 'width')

    matrix_varied_jet_norm = [normalize_cov_matrix(matrix_orig=matrix_varied_jet_orig[v], root_hist=corr_hist[4]) + matrix_norm[4] for v in range(len(var_factors))]

    corr_hist_varied_cons_norm = [corr_hist_varied_cons[v].Clone() for v in range(len(var_factors))]
    for v in range(len(var_factors)):
        corr_hist_varied_cons_norm[v].Scale(1 / corr_hist_varied_cons_norm[v].Integral(), 'width')

    matrix_varied_cons_norm = [normalize_cov_matrix(matrix_orig=matrix_varied_cons_orig[v], root_hist=corr_hist[4]) + matrix_norm[4] for v in range(len(var_factors))]

    plot_corr_hist(corr_hists=corr_hist_norm,
                   filename_graphic='chi2_plots/chi2_pt_varied_28_hist/corr_hist_norm_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                   sample_names=sample_names, title='Normalized Correlator Histograms (BW-reweighted)', hist_range=None)
    plot_corr_hist(corr_hists=[corr_hist_norm[4]]+[corr_hist_varied_jet_norm[v] for v in range(len(var_factors))],
                   filename_graphic='chi2_plots/chi2_pt_varied_28_hist/corr_hist_varied_norm_jet_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                   sample_names=['original']+var_factor_names, title='Normalized Jet-p_{T} variances',  hist_range=None)
    plot_corr_hist(corr_hists=[corr_hist_norm[4]]+[corr_hist_varied_cons_norm[v] for v in range(len(var_factors))],
                   filename_graphic='chi2_plots/chi2_pt_varied_28_hist/corr_hist_varied_norm_cons_{}_{}-{}.png'.format(level, pt_range[0], pt_range[1]),
                   sample_names=['original']+var_factor_names, title='Normalized Constituent-p_{T} variances',  hist_range=None)

    for h, sample_name in enumerate(sample_names):
        plot_matrix_in_root(matrix=matrix_norm[h], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/matrix_norm_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], sample_name),
                            hist_axis_range=hist_binning, title='Normalized Covariance Matrix Sample '+sample_name)

    for v, var_name in enumerate(var_factor_names):
        plot_matrix_in_root(matrix=matrix_varied_jet_norm[v], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/matrix_varied_jet_norm_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], var_name),
                            hist_axis_range=hist_binning, title='Normalized Covariance Matrix Jet-p_{T}-variation '+var_name)

    for v, var_name in enumerate(var_factor_names):
        plot_matrix_in_root(matrix=matrix_varied_cons_norm[v], filename_graphic='chi2_plots/chi2_pt_varied_28_matrix/matrix_varied_cons_norm_{}_{}-{}_{}.png'.format(level, pt_range[0], pt_range[1], var_name),
                            hist_axis_range=hist_binning, title='Normalized Covariance Matrix Cons-p_{T}-variation '+var_name)

    for matrix_varied_norm, var_label in [(matrix_varied_jet_norm, 'jet'), (matrix_varied_cons_norm, 'cons')]:
        chi2 = [[compute_chi2(template_hist=corr_hist_norm[h], data_hist=corr_hist_norm[4], data_cov_matrix=matrix_norm[4]) for h in range(len(sample_names))]]
        for v in range(len(var_factors)):
            chi2.append([compute_chi2(template_hist=corr_hist_norm[h], data_hist=corr_hist_norm[4], data_cov_matrix=matrix_varied_norm[v]) for h in range(len(sample_names))])

        chi2_graph = [ROOT.TGraph(len(sample_names), np.array([float(item) for item in sample_names]), np.asarray(chi2[v])) for v in range(len(chi2))]
        fit_func = ROOT.TF1('pol2_fit', 'pol2', 170, 175)
        [chi2_graph[v].Fit(fit_func, 'RQ') for v in range(len(chi2))]
        fit = [chi2_graph[v].GetFunction('pol2_fit') for v in range(len(chi2))]
        obt_top_masses = [fit[v].GetMinimumX() for v in range(len(chi2))]
        print('The calculated mass of the Top-Quark equals to {:.5f} GeV.'.format(obt_top_masses[0]))
        chi2min = [fit[v].GetMinimum() for v in range(len(chi2))]
        uncertainty_stat = abs(obt_top_masses[0] - fit[0].GetX(chi2min[0] + 1, 170, 175))
        uncertainty_tot = [abs(obt_top_masses[v] - fit[v].GetX(chi2min[v] + 1, 170, 175)) for v in range(1, len(var_factors)+1)]
        print('The total uncertainty equals {:.5f} GeV.'.format(uncertainty_tot[2]))

        uncertainties = [np.sqrt(uncertainty_tot[v]**2 - uncertainty_stat**2) for v in range(len(var_factors))]
        print('Uncertainty Up: {:.5f} GeV.'.format(uncertainties[2]))
        print('Uncertainty Down: {:.5f} GeV.'.format(uncertainties[5]))

        plot_chi2(root_graph=chi2_graph, label=['Variance: '+e for e in ['original']+var_factor_names],
                  filename='chi2_plots/chi2_variations_28/chi2_varied_{}_28_{}_{}-{}.pdf'.format(var_label, level, pt_range[0], pt_range[1]),
                  obt_top_masses=[obt_top_masses[3], obt_top_masses[6]], uncertainties=[uncertainties[2], uncertainties[5]])

        plot_uncertainties(uncertainties=uncertainties,
                           filename='chi2_plots/chi2_variations_28/chi2_varied_{}_28_uncertainties_{}_{}-{}.pdf'.format(var_label, level, pt_range[0], pt_range[1]))


"""
def main_tracker_efficiency():
    for g, level in enumerate(['Gen', 'PF']):
        uncertainties = [np.zeros((4, 3)) for _ in range(len(pt_jet_ranges))]
        chi2_graph_plot = [[[] for _ in range(3)] for _ in range(4)]
        for m, eff_deltaR in enumerate([0.01, 0.05, 0.1, 0.4]):
            for n, eff_probability in enumerate([2, 5, 10]):
                for h, sample_name in enumerate(sample_names):
                    print('Working on sample "{}" ...'.format(sample_name))
                    for k, pt_range in enumerate(pt_jet_ranges):
                        (matrices_norm[g][h][k],
                         root_hist[g][h][k],
                         root_hist_norm[g][h][k],
                         hist_range) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                            hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_{:}_{:}_{:}'.format(level, sample_name, pt_range[0], pt_range[1]),
                                                            plot_matrix=False, id_level=level, id_sample=sample_name, id_range=pt_range)

                for k, pt_jet_range in enumerate(pt_jet_ranges):
                    hists_varied[g][k][0] = prepare_histogram(filename_root_hist=filename,
                                                              hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_cons_eta_phi_{:}_None_{:}_{:}_{:}_{:}'.format(level,
                                                                                                                                                                                  pt_jet_range[0],
                                                                                                                                                                                  pt_jet_range[1],
                                                                                                                                                                                  eff_deltaR, eff_probability))

                num_bins = root_hist[g][0][0].GetNbinsX()
                for k, pt_jet_range in enumerate(pt_jet_ranges):
                    sigma[g][k][0] = [hists_varied[g][k][0].GetBinContent(i+1) - root_hist[g][4][k].GetBinContent(i+1) for i in range(num_bins)]

                    matrix_varied_orig = np.zeros((num_bins, num_bins), dtype=np.float64)
                    for i in range(num_bins):
                        for j in range(num_bins):
                            matrix_varied_orig[i, j] = sigma[g][k][0][i] * sigma[g][k][0][j]

                    matrices_varied_norm[g][k][0] = normalize_cov_matrix(matrix_orig=matrix_varied_orig, root_hist=root_hist[g][4][k]) + matrices_norm[g][4][k]

                    chi2[g][k][0] = [compute_chi2(template_hist=root_hist_norm[g][h][k], data_hist=root_hist_norm[g][4][k], data_cov_matrix=matrices_norm[g][4][k]) for h in range(9)]
                    chi2[g][k][1] = [compute_chi2(template_hist=root_hist_norm[g][h][k], data_hist=root_hist_norm[g][4][k], data_cov_matrix=matrices_varied_norm[g][k][0]) for h in range(9)]

                    chi2_graph = [ROOT.TGraph(9, np.array([171.5, 171.75, 172.0, 172.25, 172.5, 172.75, 173.0, 173.25, 173.5]), np.asarray(chi2[g][k][t])) for t in range(2)]
                    chi2_graph_plot[m][n].append(chi2_graph[1])
                    fit_func = ROOT.TF1('pol2_fit', 'pol2', 170, 175)
                    [chi2_graph[t].Fit(fit_func, 'RQ') for t in range(2)]
                    fit = [chi2_graph[t].GetFunction('pol2_fit') for t in range(2)]
                    obt_top_masses = [fit[t].GetMinimumX() for t in range(2)]
                    chi2min = [fit[t].GetMinimum() for t in range(2)]
                    uncertainty_stat = abs(obt_top_masses[0] - fit[0].GetX(chi2min[0]+1, 170, 175))
                    uncertainty_tot = abs(obt_top_masses[1] - fit[1].GetX(chi2min[1]+1, 170, 175))
                    print('The total uncertainty equals {:.5f} GeV.'.format(uncertainty_tot))

                    uncertainties[k][m, n] = np.sqrt(uncertainty_tot**2 - uncertainty_stat**2)

        uncertainties_root = [ROOT.TH2D("Uncertainties", "Systematic Uncertainties", 4, 0, 4, 3, 0, 3) for _ in range(len(pt_jet_ranges))]
        for k, pt_jet_range in enumerate(pt_jet_ranges):
            for m in range(4):
                for n in range(3):
                    uncertainties_root[k].SetBinContent(m+1, n+1, uncertainties[k][m, n])

            ROOT.gStyle.SetOptStat(0)  # Do not display stat box
            c = ROOT.TCanvas('c', 'c', 1000, 1000)
            ROOT.gPad.SetRightMargin(0.2)
            uncertainties_root[k].SetTitle('Systematic uncertainties for tracker efficiency')
            uncertainties_root[k].GetXaxis().SetTitle('Detection limit of #DeltaR')
            uncertainties_root[k].GetXaxis().SetBinLabel(1, '0.01')
            uncertainties_root[k].GetXaxis().SetBinLabel(2, '0.05')
            uncertainties_root[k].GetXaxis().SetBinLabel(3, '0.1')
            uncertainties_root[k].GetXaxis().SetBinLabel(4, '0.4')
            uncertainties_root[k].GetYaxis().SetTitle('Probability of lost events')
            uncertainties_root[k].GetYaxis().SetBinLabel(1, '50%')
            uncertainties_root[k].GetYaxis().SetBinLabel(2, '20%')
            uncertainties_root[k].GetYaxis().SetBinLabel(3, '10%')
            uncertainties_root[k].Draw('COLZ')
            c.Print(plot_directory+'chi2_plots/chi2_variations_28/chi2_tracker_efficiency_28_uncertainties_{}_{}-{}.pdf'.format(level, pt_jet_range[0], pt_jet_range[1]))

            plot_chi2(root_graph=[chi2_graph_plot[m][1][k] for m in range(3)], label=['#DeltaR: ' + e for e in ['0.01', '0.05', '0.1']],
                      filename='chi2_plots/chi2_variations_28/chi2_tracker_efficiency_deltaR_28_{}_{}-{}.pdf'.format(level, pt_jet_range[0], pt_jet_range[1]),
                      obt_top_masses=None, uncertainties=None)
            plot_chi2(root_graph=[chi2_graph_plot[1][n][k] for n in range(3)], label=['Probability: ' + e for e in ['50%', '20%', '10%']],
                      filename='chi2_plots/chi2_variations_28/chi2_tracker_efficiency_probability_28_{}_{}-{}.pdf'.format(level, pt_jet_range[0], pt_jet_range[1]),
                      obt_top_masses=None, uncertainties=None)
"""


if __name__ == '__main__':
    filename = 'histogram_files/correlator_hist_trip_28.root'
    sample_names = ['171.5', '171.75', '172.0', '172.25', '172.5', '172.75', '173.0', '173.25', '173.5']
    sample_names_stored = sample_names[:]
    sample_names_stored[sample_names_stored.index('172.5')] = 'None'

    var_factors = [1.1, 1.05, 1.02, 1.01, 0.99, 0.98, 0.95, 0.9]
    var_factor_names = ['+ 10 %', '+ 5 %', '+ 2 %', '+ 1 %', '- 1 %', '- 2 %', '- 5 %', '- 10 %']

    hist_binning = np.array([0.9, 1.2, 1.55, 1.72, 1.87, 1.9, 2.4, 2.8])

    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    for level in ['Gen']:   #  'PF']:
        for pt_range in [pt_jet_ranges[1]]:
            main(level=level, pt_range=pt_range)
            # main_tracker_efficiency(level, pt_range)

        f = ROOT.TFile(filename, 'read')
        histogram = f.Get('/Others/'+level+'-Level/number_of_events_{:}_{:}_{:}'.format(level, 450, 500))
        histogram.SetDirectory(ROOT.nullptr)
        f.Close()
        print('Number of events ('+level+'): {:.0f}'.format(histogram.GetBinContent(1)))

        f = ROOT.TFile(filename, 'read')
        histogram = f.Get('/Others/'+level+'-Level/weighted_number_of_events_{:}_{:}_{:}'.format(level, 450, 500))
        histogram.SetDirectory(ROOT.nullptr)
        f.Close()
        print('Number of events weighted by event weight ('+level+'): {:.2f}'.format(histogram.GetBinContent(1)))

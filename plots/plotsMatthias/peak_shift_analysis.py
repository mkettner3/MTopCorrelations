# -*- coding: utf-8 -*-

"""
Script to calculate the amount by which the peaks of different Top-mass-samples are shifted in relation to the peak of the data sample.
"""

import ROOT
from MTopCorrelations.Tools.user import plot_directory
from calc_triplet_and_hist_bw import rew_samples
import numpy as np


def prepare_histogram(filename_root_hist, hist_name, bin_number):
    # type: (str, str, int) -> Any

    f = ROOT.TFile(filename_root_hist, 'read')
    root_hist = f.Get(hist_name)
    root_hist.SetDirectory(ROOT.nullptr)            # Returns a pointer to root_hist in memory.
    f.Close()                                       # f.Get only returns a handle, which gets lost when TFile is closed

    group_size = 20 if bin_number is None else int(root_hist.GetNbinsX()/bin_number)
    hist_new = root_hist.Rebin(group_size, 'hist_new')

    return hist_new


def find_peak_argmax(root_hist, window):
    # type: (Any, tuple) -> tuple

    max_bin = None
    ymax = None
    for i in range(root_hist.GetNbinsX()):
        if root_hist.GetXaxis().GetBinLowEdge(i+1) > window[0] and root_hist.GetXaxis().GetBinUpEdge(i+1) < window[1]:
            bin_content = root_hist.GetBinContent(i+1)
            if bin_content > ymax:
                ymax = bin_content
                max_bin = i+1
    max_value = root_hist.GetBinCenter(max_bin) if max_bin is not None else None

    return max_bin, max_value


def perform_gauss_fit(root_hist, window):
    # type: (Any, tuple) -> tuple

    gaussFitFunction = ROOT.TF1("gaussFitFunction", "gaus", window[0], window[1])
    root_hist.Fit(gaussFitFunction, "RQ")
    gauss_fit = root_hist.GetFunction("gaussFitFunction")

    norm_gauss = gauss_fit.GetParameter(0)
    mean_gauss = gauss_fit.GetParameter(1)
    sigma_gauss = gauss_fit.GetParameter(2)
    max_gauss = gauss_fit.GetMaximumX()

    return norm_gauss, mean_gauss, sigma_gauss, max_gauss


def perform_general_crystal_ball_fit(root_hist, window, norm_gauss, mean_gauss, sigma_gauss, double_sided=False):
    # type: (Any, tuple, float, float, float, bool) -> tuple

    norm_init = norm_gauss
    mean_init = mean_gauss
    sigma_init = sigma_gauss
    alphaL_init = 2 * sigma_gauss
    alphaR_init = 2 * sigma_gauss
    nL_init = 2
    nR_init = 2

    if double_sided:
        CBFunction = ROOT.TF1("CBFunction", DoubleSidedCB, window[0], window[1], 8)
        CBFunction.SetParameters(norm_init, mean_init, sigma_init, sigma_init, alphaL_init, alphaR_init, nL_init, nR_init)
        CBFunction.SetParLimits(2, 0.3*sigma_gauss, 10)
        CBFunction.SetParLimits(3, 0.3*sigma_gauss, 10)
        CBFunction.SetParLimits(6, 0, 20)
        CBFunction.SetParLimits(7, 0, 20)
    else:
        CBFunction = ROOT.TF1("CBFunction", OneSidedCB, window[0], window[1], 7)
        CBFunction.SetParameters(norm_init, mean_init, sigma_init, alphaL_init, alphaR_init, nL_init, nR_init)
        CBFunction.SetParLimits(2, 0.3*sigma_gauss, 10)
        CBFunction.SetParLimits(5, 0, 20)
        CBFunction.SetParLimits(6, 0, 20)

    root_hist.Fit(CBFunction, "RQ")
    CB_fit = root_hist.GetFunction("CBFunction")

    crystal_ball_mean = CB_fit.GetParameter(1)
    crystal_ball_mean_error = CB_fit.GetParError(1)
    crystal_ball_max = CB_fit.GetMaximumX()
    standard_parameters = [CB_fit.GetParameter(i) for i in range(2, CB_fit.GetNpar())]

    print(CB_fit.GetParameter(3), CB_fit.GetParameter(4))

    return crystal_ball_mean, crystal_ball_mean_error, crystal_ball_max, standard_parameters


def perform_standardized_crystal_ball_fit(root_hist, window, norm_gauss, mean_gauss):
    # type: (Any, tuple, float, float) -> tuple

    CBFunction = ROOT.TF1("CBFunction", OneSidedCB_standard, window[0], window[1], 2)
    CBFunction.SetParameters(norm_gauss, mean_gauss)

    root_hist.Fit(CBFunction, "RQ")
    CB_fit = root_hist.GetFunction("CBFunction")

    crystal_ball_mean = CB_fit.GetParameter(1)
    crystal_ball_mean_error = CB_fit.GetParError(1)
    crystal_ball_max = CB_fit.GetMaximumX()
    ndf = CB_fit.GetNDF()
    chi2 = CB_fit.GetChisquare()

    return crystal_ball_mean, crystal_ball_mean_error, crystal_ball_max, chi2, chi2/ndf


def DoubleSidedCB(x, params):
    norm, mu, widthL, widthR, alphaL, alphaR, nL, nR = params

    signAL = 1.0 if nL > 0.0 else -1.
    signAR = 1.0 if nR > 0.0 else -1.

    AL = signAL * pow(abs(nL)/abs(alphaL), nL) * np.exp(-abs(alphaL**2)/2.0)
    AR = signAR * pow(abs(nR)/abs(alphaR), nR) * np.exp(-abs(alphaR**2)/2.0)
    BL = nL/abs(alphaL) - abs(alphaL)
    BR = nR/abs(alphaR) - abs(alphaR)

    diffL = (x[0]-mu)/widthL
    diffR = (x[0]-mu)/widthR

    if diffL <= -alphaL:
        result = AL * pow(abs(BL-diffL), -nL)
    elif -alphaL < diffL < 0.:
        result = np.exp(-0.5 * diffL**2)
    elif 0. < diffR < alphaR:
        result = np.exp(-0.5 * diffR**2)
    elif diffR >= alphaR:
        result = AR * pow(abs(BR+diffR), -nR)
    else:
        raise RuntimeError('Crystal-Ball-function is out of bounds! ({}, {}, {}, {})'.format(diffL, diffR, alphaL, alphaR))

    return norm*result


def OneSidedCB(x, params):
    norm, mu, width, alphaL, alphaR, nL, nR = params

    AL = pow(nL/abs(alphaL), nL) * np.exp(-abs(alphaL**2)/2.0)
    AR = pow(nR/abs(alphaR), nR) * np.exp(-abs(alphaR**2)/2.0)
    BL = nL/abs(alphaL) - abs(alphaL)
    BR = nR/abs(alphaR) - abs(alphaR)

    diff = (x[0]-mu)/width

    if diff <= -alphaL:
        result = AL * pow(abs(BL-diff), -nL)
    elif -alphaL < diff < alphaR:
        result = np.exp(-0.5 * diff**2)
    elif diff >= alphaR:
        result = AR * pow(abs(BR+diff), -nR)
    else:
        raise RuntimeError('Crystal-Ball-function is out of bounds! ({}, {}, {})'.format(diff, alphaL, alphaR))

    return norm*result


def OneSidedCB_standard(x, params):
    norm, mu = params
    width, alphaL, alphaR, nL, nR = [0.15152, 0.38339, 1.01128, 20, 20]

    AL = pow(nL/abs(alphaL), nL) * np.exp(-abs(alphaL**2)/2.0)
    AR = pow(nR/abs(alphaR), nR) * np.exp(-abs(alphaR**2)/2.0)
    BL = nL/abs(alphaL) - abs(alphaL)
    BR = nR/abs(alphaR) - abs(alphaR)

    diff = (x[0]-mu)/width

    if diff <= -alphaL:
        result = AL * pow(abs(BL-diff), -nL)
    elif -alphaL < diff < alphaR:
        result = np.exp(-0.5 * diff**2)
    elif diff >= alphaR:
        result = AR * pow(abs(BR+diff), -nR)
    else:
        raise RuntimeError('Crystal-Ball-function is out of bounds! ({}, {}, {})'.format(diff, alphaL, alphaR))

    return norm*result


def draw_histogram(root_hist, filename_graphic, sample_name, fit_label, fit_function, chi2_ndf=None, ylim=(0, 0.00015)):
    # type: (Any, str, str, str, str, float, tuple) -> None

    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)  # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)  # Axis ticks right
    ROOT.gStyle.SetOptStat(0)  # Do not display stat box

    c = ROOT.TCanvas('c', 'c', 600, 600)
    legend = ROOT.TLegend(0.65, 0.77, 0.85, 0.87)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)

    root_hist.SetLineColor(ROOT.kBlack)
    root_hist.SetTitle('')
    root_hist.SetLineWidth(2)
    root_hist.SetLineStyle(1)
    legend.AddEntry(root_hist, 'm_{top} = '+sample_name+' GeV', 'l')
    root_hist.GetFunction(fit_function).SetLineColor(ROOT.kRed)
    legend.AddEntry(root_hist.GetFunction(fit_function), fit_label, 'l')
    if chi2_ndf is not None:
        legend.AddEntry(ROOT.nullptr, '#chi^{2}/NDF = '+'{:.2f}'.format(chi2_ndf), '')

    root_hist.GetXaxis().SetRangeUser(0, 2.5)  # x-axis range (also works for y-axis)
    root_hist.GetXaxis().SetTitle("3#zeta")
    root_hist.GetXaxis().CenterTitle(ROOT.kTRUE)
    root_hist.GetXaxis().SetNdivisions(5, 5, 0)  # Unterteilung der x-Achse
    root_hist.GetXaxis().SetTitleOffset(1.5)
    root_hist.GetYaxis().SetRangeUser(ylim[0], ylim[1])
    # root_hist.GetYaxis().SetTitle("Energy-weighted Triplets")
    root_hist.GetYaxis().CenterTitle(ROOT.kTRUE)
    root_hist.GetYaxis().SetNdivisions(5, 5, 0)    # 2*1000000 + (root_hist.GetYaxis().GetNdiv() % 1000000))  # Unterteilung der y-Achse
    root_hist.GetYaxis().SetMaxDigits(3)  # 3 ist die einzig sinnvolle Einstellung, weil als Exponent der Zehnerpotenz nur Vielfache von 3 verwendet werden.
    root_hist.GetYaxis().SetTitleOffset(1.5)

    root_hist.Draw('HIST')
    root_hist.GetFunction(fit_function).Draw('SAME')
    legend.Draw()
    c.Print(plot_directory+filename_graphic)


def draw_mass_fit_graph(correlator_values_variations, filename_graphic, chart_title, correlator_value_errors=None):
    # type: (list, str, str, list) -> dict

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    legend = ROOT.TLegend(0.15, 0.62, 0.28, 0.87)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)

    top_mass_graph = {}
    for s, (label, correlator_values) in enumerate(correlator_values_variations):
        if correlator_value_errors is None:
            top_mass_graph[s] = ROOT.TGraph(len(correlator_values), np.array(list(correlator_values.keys())), np.asarray(list(correlator_values.values())))
        else:
            top_mass_graph[s] = ROOT.TGraphErrors(len(correlator_values), np.array(list(correlator_values.keys())), np.asarray(list(correlator_values.values())),
                                                  np.zeros(len(correlator_values)), np.array(list(correlator_value_errors[s].values())))
        fit_func = ROOT.TF1('pol2_fit', 'pol2', 150, 195)
        top_mass_graph[s].Fit(fit_func, 'RQ', '', 170, 175)

        top_mass_graph[s].SetMarkerSize(2)
        top_mass_graph[s].SetMarkerStyle(47)
        top_mass_graph[s].SetMarkerColor(1)
        top_mass_graph[s].GetFunction('pol2_fit').SetLineColor(s+1)
        legend.AddEntry(top_mass_graph[s].GetFunction('pol2_fit'), label, 'l')

    top_mass_graph[0].SetTitle(chart_title)
    top_mass_graph[0].GetXaxis().SetTitle('Top-Mass (GeV)')
    top_mass_graph[0].GetXaxis().CenterTitle(ROOT.kTRUE)
    top_mass_graph[0].GetXaxis().SetNdivisions(505)  # Unterteilung der x-Achse
    top_mass_graph[0].GetXaxis().SetTitleOffset(1.5)
    top_mass_graph[0].GetYaxis().SetTitle("3#zeta")
    top_mass_graph[0].GetYaxis().CenterTitle(ROOT.kTRUE)
    top_mass_graph[0].GetYaxis().SetNdivisions(505)
    top_mass_graph[0].GetYaxis().SetTitleOffset(1.5)
    top_mass_graph[0].GetYaxis().SetRangeUser(min([top_mass_graph[i].GetHistogram().GetMinimum() for i in top_mass_graph]),
                                              max([top_mass_graph[i].GetHistogram().GetMaximum() for i in top_mass_graph]))
    top_mass_graph[0].Draw('AP')
    for s in range(len(correlator_values_variations)):
        top_mass_graph[s].Draw('P SAME')        # This object must not be deleted until c.Print() is executed.
    legend.Draw()
    c.Print(plot_directory+filename_graphic)

    return top_mass_graph


def generate_fits(rew_samples, mtop_bw_names, filename, subfolder, pt_jet_range, peak_window, variation_id='', bin_number=None):
    # type: (list, list, str, str, tuple, tuple, str, int) -> tuple

    # hist = prepare_histogram(filename_root_hist=filename,
    #                          hist_name='/Top-Quark/Gen-Level/weighted/correlator_hist_Gen_{:}_{:}_{:}'.format(None, pt_jet_range[0], pt_jet_range[1]),
    #                          bin_number=bin_number)
    # gauss_norm, gauss_mean, gauss_sigma, _ = perform_gauss_fit(root_hist=hist, window=peak_window)
    # _, _, _, standard_param = perform_general_crystal_ball_fit(root_hist=hist, window=peak_window, norm_gauss=gauss_norm,
    #                                                            mean_gauss=gauss_mean, sigma_gauss=gauss_sigma, double_sided=False)
    # print(standard_param)

    gauss_means = {}
    gauss_maximums = {}
    CB_means = {}
    CB_mean_errors = {}
    CB_maximums = {}

    for mtop_bw, mtop_bw_name in zip(rew_samples, mtop_bw_names):
        hist = prepare_histogram(filename_root_hist=filename,
                                 hist_name='/Top-Quark/Gen-Level/weighted/correlator_hist{:}_Gen_{:}_{:}_{:}'.format(variation_id, mtop_bw, pt_jet_range[0], pt_jet_range[1]),
                                 bin_number=bin_number)

        # Find the maximum bin
        peak_max_bin, peak_max_value = find_peak_argmax(root_hist=hist, window=peak_window)
        print('Maximum bin number: {:}; Maximum value: {:}'.format(peak_max_bin, peak_max_value))

        # Perform a Gauss fit to get a rough idea of fit parameters
        gauss_norm, gauss_mean, gauss_sigma, gauss_max = perform_gauss_fit(root_hist=hist, window=peak_window)
        print('Gauss-Fit mean: {:.3f}; Gauss-Fit Sigma: {:.3f}'.format(gauss_mean, gauss_sigma))
        draw_histogram(root_hist=hist, filename_graphic=subfolder+'/peak_fitting/correlator_gauss_fit{:}_Gen_{:}_{:}-{:}.pdf'.format(variation_id, mtop_bw_name, pt_jet_range[0], pt_jet_range[1]),
                       sample_name=mtop_bw_name, fit_label='Gauss-Fit', fit_function='gaussFitFunction')

        # Now fit with double-sided crystal ball function
        (CB_mean, CB_mean_error, CB_max,
         CB_chi2, CB_chi2_ndf) = perform_standardized_crystal_ball_fit(root_hist=hist, window=peak_window, norm_gauss=gauss_norm, mean_gauss=gauss_mean)
        print('Crystal-Ball-Fit mean: {:.3f} +- {:.3f}; Chi2: {:.3f}; Chi2/NDF: {:.3f}'.format(CB_mean, CB_mean_error, CB_chi2, CB_chi2_ndf))
        draw_histogram(root_hist=hist, filename_graphic=subfolder+'/peak_fitting/correlator_CB_fit{:}_Gen_{:}_{:}-{:}.pdf'.format(variation_id, mtop_bw_name, pt_jet_range[0], pt_jet_range[1]),
                       sample_name=mtop_bw_name, fit_label='Crystal-Ball-Fit', fit_function='CBFunction', chi2_ndf=CB_chi2_ndf)

        gauss_means[mtop_bw if mtop_bw is not None else 172.5] = gauss_mean
        gauss_maximums[mtop_bw if mtop_bw is not None else 172.5] = gauss_max
        CB_means[mtop_bw if mtop_bw is not None else 172.5] = CB_mean
        CB_mean_errors[mtop_bw if mtop_bw is not None else 172.5] = CB_mean_error
        CB_maximums[mtop_bw if mtop_bw is not None else 172.5] = CB_max

        print(' ')

    return gauss_means, gauss_maximums, CB_means, CB_mean_errors, CB_maximums


def calc_variation_shifts(mass_fit_graph, var_values, var_corr_values, filename_graphic):
    # type: (Any, list, list, str) -> dict

    std_val = (var_values[0]+var_values[-1])/2

    var_shifts = {}
    for k, var_value in enumerate([std_val]+var_values):
        var_shifts[var_value] = mass_fit_graph.GetFunction('pol2_fit').GetX(var_corr_values[k][1][172.5])

    # print('====    ====    ====    ====')
    # print(original_corr_value)
    # for k in [1]+var_values:
    #     print(k, var_shifts[k])
    var_shifts = {var_value: var_shifts[var_value] - 172.5 for var_value in [std_val]+var_values}

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)

    var_shifts_graph = ROOT.TGraph(len(var_shifts), np.asarray(list(var_shifts.keys())), np.asarray(list(var_shifts.values())))

    var_shifts_graph.SetMarkerSize(2)
    var_shifts_graph.SetMarkerStyle(47)
    var_shifts_graph.SetMarkerColor(1)
    var_shifts_graph.SetTitle('Variation Shifts')
    var_shifts_graph.GetXaxis().SetTitle('Variation Factor')
    var_shifts_graph.GetXaxis().CenterTitle(ROOT.kTRUE)
    var_shifts_graph.GetXaxis().SetNdivisions(505)  # Unterteilung der x-Achse
    var_shifts_graph.GetXaxis().SetTitleOffset(1.5)
    var_shifts_graph.GetYaxis().SetTitle('Shift in the Top Mass (GeV)')
    var_shifts_graph.GetYaxis().CenterTitle(ROOT.kTRUE)
    var_shifts_graph.GetYaxis().SetNdivisions(505)
    var_shifts_graph.GetYaxis().SetTitleOffset(1.5)
    var_shifts_graph.Draw('AP')

    c.Print(plot_directory+filename_graphic)

    return var_shifts


def calc_efficiency_shifts(mass_fit_graph, var_values, var_corr_values, filename_graphic):
    # type: (Any, list, list, str) -> dict

    var_shifts = {}
    k = 0
    for m in var_values[0]:
        for n in var_values[1]:
            var_shifts[(m, n)] = mass_fit_graph.GetFunction('pol2_fit').GetX(var_corr_values[k+1][1][172.5])
            k += 1

    # print('====    ====    ====    ====')
    # print(original_corr_value)
    # for k in [1]+var_values:
    #     print(k, var_shifts[k])
    var_shifts = {(m, n): var_shifts[(m, n)] - 172.5 for m in var_values[0] for n in var_values[1]}

    efficiency_graph = ROOT.TH2D("Efficiencies", "Efficiency Variations", len(var_values[0]), 0, len(var_values[0]), len(var_values[1]), 0, len(var_values[1]))
    for i, m in enumerate(var_values[0]):
        for j, n in enumerate(var_values[1]):
            efficiency_graph.SetBinContent(i+1, j+1, var_shifts[(m, n)])

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetRightMargin(0.2)
    efficiency_graph.SetTitle('Systematic uncertainties for tracker efficiency')
    efficiency_graph.GetXaxis().SetTitle('Detection limit of #DeltaR')
    for i, m in enumerate(var_values[0]):
        efficiency_graph.GetXaxis().SetBinLabel(i+1, '{:}'.format(m))
    efficiency_graph.GetYaxis().SetTitle('Probability of lost events')
    for j, n in enumerate(var_values[1]):
        efficiency_graph.GetYaxis().SetBinLabel(j+1, '{:}%'.format(n))
    efficiency_graph.Draw('COLZ')
    c.Print(plot_directory+filename_graphic)

    return var_shifts


if __name__ == '__main__':
    subfolder = '/generation_33'
    filename = 'histogram_files/correlator_hist_trip_33.root'
    pt_jet_range = (450, 500)
    peak_window = (0.7, 1.7)
    mtop_bw_names = ['{:.2f}'.format(elem) if elem is not None else '172.50' for elem in rew_samples]
    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    jet_pt_variations = [1.02, 1.01, 1.005, 1.002, 1.001, 0.999, 0.998, 0.995, 0.99, 0.98]
    jet_pt_var_ids = ['_varied_jet_{:.3f}'.format(v) for v in jet_pt_variations]
    jet_pt_var_names = ['+ 2 %', '+ 1 %', '+ 0.5 %', '+ 0.2 %', '+ 0.1 %', '- 0.1 %', '- 0.2 %', '- 0.5 %', '- 1 %', '- 2 %']

    cons_pt_variations = [-2, -1, -0.5, 0.5, 1, 2]
    cons_pt_var_ids = ['_varied_cons_pt_{:.2f}'.format(v) for v in cons_pt_variations]
    cons_pt_var_names = ['{:.1f}'.format(v) for v in cons_pt_variations]

    deltaR_variations = [0.01, 0.05, 0.1]
    probab_variations = [2, 5, 10]
    efficiency_variations = [deltaR_variations, [50, 20, 10]]
    efficiency_ids = ['_varied_efficiency_{:.2f}_{:.2f}'.format(eff_deltaR, eff_probability) for eff_deltaR in deltaR_variations for eff_probability in probab_variations]
    efficiency_names = ['{:}, {:}'.format(delR, probab) for delR in deltaR_variations for probab in ['50 %', '20 %', '10 %']]

    # ========== Variations of Jet-pT ========== #
    gauss_means_all_jet_pt_var = []
    gauss_maximums_all_jet_pt_var = []
    CB_means_all_jet_pt_var = []
    CB_maximums_all_jet_pt_var = []
    for var_id, var_name in zip(['']+jet_pt_var_ids, ['original']+jet_pt_var_names):
        (gauss_means, gauss_maximums,
         CB_means, CB_mean_errors, CB_maximums) = generate_fits(rew_samples=rew_samples, mtop_bw_names=mtop_bw_names, filename=filename, subfolder=subfolder,
                                                                pt_jet_range=pt_jet_range, peak_window=peak_window, variation_id=var_id)
        gauss_means_all_jet_pt_var.append((var_name, gauss_means))
        gauss_maximums_all_jet_pt_var.append((var_name, gauss_maximums))
        CB_means_all_jet_pt_var.append((var_name, CB_means))
        CB_maximums_all_jet_pt_var.append((var_name, CB_maximums))

    draw_mass_fit_graph(correlator_values_variations=gauss_means_all_jet_pt_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_fit_jet.pdf', chart_title='Top Mass Gauss Fit (Means)')
    draw_mass_fit_graph(correlator_values_variations=gauss_maximums_all_jet_pt_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_max_fit_jet.pdf', chart_title='Top Mass Gauss Fit (Maximums)')
    fit_graph_CB_means = draw_mass_fit_graph(correlator_values_variations=CB_means_all_jet_pt_var,
                                             filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit_jet.pdf', chart_title='Top Mass Crystal Ball Fit (Means)') # , correlator_value_errors=CB_mean_errors
    draw_mass_fit_graph(correlator_values_variations=CB_maximums_all_jet_pt_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_CB_max_fit_jet.pdf', chart_title='Top Mass Crystal Ball Fit (Maximums)')

    jet_pt_var_shifts = calc_variation_shifts(mass_fit_graph=fit_graph_CB_means[0], var_values=jet_pt_variations,
                                              var_corr_values=CB_means_all_jet_pt_var, filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit_jet_shifts.pdf')

    # ========== Variations of Constituent-pT ========== #
    gauss_means_all_cons_pt_var = []
    gauss_maximums_all_cons_pt_var = []
    CB_means_all_cons_pt_var = []
    CB_maximums_all_cons_pt_var = []
    for var_id, var_name in zip(['']+cons_pt_var_ids, ['original']+cons_pt_var_names):
        (gauss_means, gauss_maximums,
         CB_means, CB_mean_errors, CB_maximums) = generate_fits(rew_samples=rew_samples, mtop_bw_names=mtop_bw_names, filename=filename, subfolder=subfolder,
                                                                pt_jet_range=pt_jet_range, peak_window=peak_window, variation_id=var_id)
        gauss_means_all_cons_pt_var.append((var_name, gauss_means))
        gauss_maximums_all_cons_pt_var.append((var_name, gauss_maximums))
        CB_means_all_cons_pt_var.append((var_name, CB_means))
        CB_maximums_all_cons_pt_var.append((var_name, CB_maximums))

    draw_mass_fit_graph(correlator_values_variations=gauss_means_all_cons_pt_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_fit_cons.pdf', chart_title='Top Mass Gauss Fit (Means)')
    draw_mass_fit_graph(correlator_values_variations=gauss_maximums_all_cons_pt_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_max_fit_cons.pdf', chart_title='Top Mass Gauss Fit (Maximums)')
    fit_graph_CB_means = draw_mass_fit_graph(correlator_values_variations=CB_means_all_cons_pt_var,
                                             filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit_cons.pdf', chart_title='Top Mass Crystal Ball Fit (Means)') # , correlator_value_errors=CB_mean_errors
    draw_mass_fit_graph(correlator_values_variations=CB_maximums_all_cons_pt_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_CB_max_fit_cons.pdf', chart_title='Top Mass Crystal Ball Fit (Maximums)')

    jet_cons_var_shifts = calc_variation_shifts(mass_fit_graph=fit_graph_CB_means[0], var_values=cons_pt_variations,
                                                var_corr_values=CB_means_all_cons_pt_var, filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit_cons_shifts.pdf')

    # ========== Efficiency analysis ========== #
    gauss_means_all_efficiency_var = []
    gauss_maximums_all_efficiency_var = []
    CB_means_all_efficiency_var = []
    CB_maximums_all_efficiency_var = []
    for var_id, var_name in zip(['']+efficiency_ids, ['original']+efficiency_names):
        (gauss_means, gauss_maximums,
         CB_means, CB_mean_errors, CB_maximums) = generate_fits(rew_samples=rew_samples, mtop_bw_names=mtop_bw_names, filename=filename, subfolder=subfolder,
                                                                pt_jet_range=pt_jet_range, peak_window=peak_window, variation_id=var_id)
        gauss_means_all_efficiency_var.append((var_name, gauss_means))
        gauss_maximums_all_efficiency_var.append((var_name, gauss_maximums))
        CB_means_all_efficiency_var.append((var_name, CB_means))
        CB_maximums_all_efficiency_var.append((var_name, CB_maximums))

    draw_mass_fit_graph(correlator_values_variations=gauss_means_all_efficiency_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_fit_efficiency.pdf', chart_title='Top Mass Gauss Fit (Means)')
    draw_mass_fit_graph(correlator_values_variations=gauss_maximums_all_efficiency_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_max_fit_efficiency.pdf', chart_title='Top Mass Gauss Fit (Maximums)')
    fit_graph_CB_means = draw_mass_fit_graph(correlator_values_variations=CB_means_all_efficiency_var,
                                             filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit_efficiency.pdf',
                                             chart_title='Top Mass Crystal Ball Fit (Means)')  # , correlator_value_errors=CB_mean_errors
    draw_mass_fit_graph(correlator_values_variations=CB_maximums_all_efficiency_var,
                        filename_graphic=subfolder+'/peak_fitting/top_mass_CB_max_fit_efficiency.pdf', chart_title='Top Mass Crystal Ball Fit (Maximums)')

    efficiency_var_shifts = calc_efficiency_shifts(mass_fit_graph=fit_graph_CB_means[0], var_values=efficiency_variations,
                                                   var_corr_values=CB_means_all_efficiency_var, filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit_efficiency_shifts.pdf')

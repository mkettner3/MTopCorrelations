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

    group_size = int(root_hist.GetNbinsX()/bin_number)
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

    return gauss_fit, norm_gauss, mean_gauss, sigma_gauss


def perform_crystal_ball_fit(root_hist, window, norm_gauss, mean_gauss, sigma_gauss, double_sided=True):
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
        CBFunction.SetParLimits(2, 0, 10)
        CBFunction.SetParLimits(3, 0, 10)
        CBFunction.SetParLimits(6, 0, 5)
        CBFunction.SetParLimits(7, 0, 5)
    else:
        CBFunction = ROOT.TF1("CBFunction", OneSidedCB, window[0], window[1], 7)
        CBFunction.SetParameters(norm_init, mean_init, sigma_init, alphaL_init, alphaR_init, nL_init, nR_init)
        CBFunction.SetParLimits(2, 0, 10)
        CBFunction.SetParLimits(5, 0, 5)
        CBFunction.SetParLimits(6, 0, 5)

    root_hist.Fit(CBFunction, "RQ")
    CB_fit = root_hist.GetFunction("CBFunction")

    crystal_ball_mean = CB_fit.GetParameter(1)

    return CB_fit, crystal_ball_mean


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


def draw_histogram(root_hist, filename_graphic, sample_name, fit_label, fit_function, ylim=(0, 0.0001)):
    # type: (Any, str, str, str, str, tuple) -> None

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

    root_hist.GetXaxis().SetRangeUser(0, 2.5)  # x-axis range (also works for y-axis)
    root_hist.GetXaxis().SetTitle("3#zeta")
    root_hist.GetXaxis().CenterTitle(ROOT.kTRUE)
    root_hist.GetXaxis().SetNdivisions(505)  # Unterteilung der x-Achse
    root_hist.GetXaxis().SetTitleOffset(1.5)
    root_hist.GetYaxis().SetRangeUser(ylim[0], ylim[1])
    root_hist.GetYaxis().SetTitle("Energy-weighted Triplets")
    root_hist.GetYaxis().CenterTitle(ROOT.kTRUE)
    root_hist.GetYaxis().SetNdivisions(505)    # 2*1000000 + (root_hist.GetYaxis().GetNdiv() % 1000000))  # Unterteilung der y-Achse
    root_hist.GetYaxis().SetMaxDigits(3)  # 3 ist die einzig sinnvolle Einstellung, weil als Exponent der Zehnerpotenz nur Vielfache von 3 verwendet werden.
    root_hist.GetYaxis().SetTitleOffset(1.5)

    root_hist.Draw('HIST')
    root_hist.GetFunction(fit_function).Draw('SAME')
    legend.Draw()
    c.Print(plot_directory+filename_graphic)


def draw_mass_fit_graph(correlator_values, filename_graphic, chart_title):
    top_mass_graph = ROOT.TGraph(len(correlator_values), np.array(list(correlator_values.keys())), np.asarray(list(correlator_values.values())))
    fit_func = ROOT.TF1('pol2_fit', 'pol2', 170, 175)
    top_mass_graph.Fit(fit_func, 'RQ')

    ROOT.gStyle.SetOptStat(0)  # Do not display stat box
    ROOT.gStyle.SetLegendBorderSize(1)  # No border for legend
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(0)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    # legend = ROOT.TLegend(0.65, 0.5, 0.85, 0.87)
    ROOT.gPad.SetLeftMargin(0.12)
    ROOT.gPad.SetBottomMargin(0.12)

    top_mass_graph.SetMarkerSize(2)
    top_mass_graph.SetMarkerStyle(47)
    top_mass_graph.SetMarkerColor(1)
    top_mass_graph.GetFunction('pol2_fit').SetLineColor(ROOT.kRed)
    # legend.AddEntry(top_mass_graph.GetFunction('pol2_fit'), label[s], 'l')
    top_mass_graph.SetTitle(chart_title)
    top_mass_graph.GetXaxis().SetTitle('Top-Mass (GeV)')
    top_mass_graph.GetXaxis().CenterTitle(ROOT.kTRUE)
    top_mass_graph.GetXaxis().SetNdivisions(505)  # Unterteilung der x-Achse
    top_mass_graph.GetXaxis().SetTitleOffset(1.5)
    top_mass_graph.GetYaxis().SetTitle("3#zeta")
    top_mass_graph.GetYaxis().CenterTitle(ROOT.kTRUE)
    top_mass_graph.GetYaxis().SetNdivisions(505)
    top_mass_graph.GetYaxis().SetTitleOffset(1.5)
    # top_mass_graph.SetMaximum(top_mass_graph.GetHistogram().GetMaximum())
    # top_mass_graph.SetMinimum(top_mass_graph.GetHistogram().GetMinimum())
    top_mass_graph.Draw('AP')
    # for g in root_graph:
    #     g.Draw('P SAME')
    # legend.Draw()

    c.Print(plot_directory+filename_graphic)


if __name__ == '__main__':
    subfolder = '/generation_32'
    filename = 'histogram_files/correlator_hist_trip_32.root'
    pt_jet_range = (450, 500)
    peak_window = (0.7, 1.7)

    mtop_bw_names = ['{:.2f}'.format(elem) if elem is not None else '172.50' for elem in rew_samples]
    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    gauss_means = {}
    CB_means = {}
    for mtop_bw, mtop_bw_name in zip(rew_samples[1:-1], mtop_bw_names[1:-1]):
        hist = prepare_histogram(filename_root_hist=filename,
                                 hist_name='/Top-Quark/Gen-Level/weighted/correlator_hist_Gen_{:}_{:}_{:}'.format(mtop_bw, pt_jet_range[0], pt_jet_range[1]),
                                 bin_number=45)

        # Find the maximum bin
        peak_max_bin, peak_max_value = find_peak_argmax(root_hist=hist, window=peak_window)
        print('Maximum bin number: {:}; Maximum value: {:}'.format(peak_max_bin,  peak_max_value))

        # Perform a Gauss fit to get a rough idea of fit parameters
        gauss_fit_func, gauss_norm, gauss_mean, gauss_sigma = perform_gauss_fit(root_hist=hist, window=peak_window)
        print('Gauss-Fit mean: {:.3f}; Gauss-Fit Sigma: {:.3f}'.format(gauss_mean, gauss_sigma))
        draw_histogram(root_hist=hist, filename_graphic=subfolder+'/peak_fitting/correlator_gauss_fit_Gen_{:}_{:}-{:}.pdf'.format(mtop_bw_name, pt_jet_range[0], pt_jet_range[1]),
                       sample_name=mtop_bw_name, fit_label='Gauss-Fit', fit_function='gaussFitFunction')

        # Now fit with double-sided crystal ball function
        CB_fit_func, CB_mean = perform_crystal_ball_fit(root_hist=hist, window=peak_window, norm_gauss=gauss_norm,
                                                        mean_gauss=gauss_mean, sigma_gauss=gauss_sigma, double_sided=False)
        print('Crystal-Ball-Fit mean: {:.3f}'.format(CB_mean))
        draw_histogram(root_hist=hist, filename_graphic=subfolder+'/peak_fitting/correlator_CB_fit_Gen_{:}_{:}-{:}.pdf'.format(mtop_bw_name, pt_jet_range[0], pt_jet_range[1]),
                       sample_name=mtop_bw_name, fit_label='Crystal-Ball-Fit', fit_function='CBFunction')
        gauss_means[mtop_bw if mtop_bw is not None else 172.5] = gauss_mean
        CB_means[mtop_bw if mtop_bw is not None else 172.5] = CB_mean

        print(' ')

    print(CB_means)
    draw_mass_fit_graph(correlator_values=gauss_means, filename_graphic=subfolder+'/peak_fitting/top_mass_gauss_fit.pdf', chart_title='Top Mass Gauss Fit')
    draw_mass_fit_graph(correlator_values=CB_means, filename_graphic=subfolder+'/peak_fitting/top_mass_CB_fit.pdf', chart_title='Top Mass Crystal Ball Fit')

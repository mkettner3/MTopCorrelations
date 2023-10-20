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


def perform_crystal_ball_fit(root_hist, window, norm_gauss, mean_gauss, sigma_gauss):
    # type: (Any, tuple, float, float, float) -> tuple

    norm_init = norm_gauss
    mean_init = mean_gauss
    sigmaL_init = sigma_gauss
    sigmaR_init = sigma_gauss
    alphaL_init = 2 * sigma_gauss
    alphaR_init = 2 * sigma_gauss
    nL_init = 2
    nR_init = 2
    rangeGaus_lo = window[0]
    rangeGaus_hi = window[1]

    CBFunction = ROOT.TF1("CBFunction", DoubleSidedCB, rangeGaus_lo, rangeGaus_hi, 8)
    CBFunction.SetParameters(norm_init, mean_init, sigmaL_init, sigmaR_init, alphaL_init, alphaR_init, nL_init, nR_init)
    CBFunction.SetParLimits(6, -5, 5)
    CBFunction.SetParLimits(7, -5, 5)

    root_hist.Fit(CBFunction, "RQ")
    CB_fit = root_hist.GetFunction("CBFunction")

    return CB_fit


def DoubleSidedCB(x, params):
    norm = params[0]
    mu = params[1]
    widthL = abs(params[2])
    widthR = abs(params[3])
    alphaL = params[4]
    alphaR = params[5]
    nL = params[6]
    nR = params[7]

    signAL = 1.0 if nL > 0.0 else -1.
    signAR = 1.0 if nR > 0.0 else -1.

    # print mu, widthL, widthR, alphaL, alphaR, nL, nR
    AL = signAL * pow(abs(nL)/abs(alphaL), nL) * np.exp(-abs(alphaL*alphaL)/2.0)
    AR = signAR * pow(abs(nR)/abs(alphaR), nR) * np.exp(-abs(alphaR*alphaR)/2.0)
    BL = nL/abs(alphaL) - abs(alphaL)
    BR = nR/abs(alphaR) - abs(alphaR)

    diffL = (x[0]-mu)/widthL
    diffR = (x[0]-mu)/widthR

    if diffL < -alphaL:
        result = AL * pow(abs(BL-diffL), -nL)
    elif -alphaL < diffL < 0.:
        result = np.exp(-0.5 * diffL*diffL)
    elif 0. < diffR < alphaR:
        result = np.exp(-0.5 * diffR*diffR)
    elif diffR > alphaR:
        result = AR * pow(abs(BR+diffR), -nR)
    else:
        raise RuntimeError('Crystal-Ball-function is out of bounds! ({}, {}, {}, {})'.format(diffL, diffR, alphaL, alphaR))

    return norm*result


def OneSidedCB(x, params):
    norm = params[0]
    mu = params[1]
    width = abs(params[2])
    alphaL = params[3]
    alphaR = params[4]
    nL = params[5]
    nR = params[6]

    signAL = 1.0
    if nL < 0.0 and not (nL % 2) == 0:
        signAL = -1

    signAR = 1.0
    if nR < 0.0 and not (nR % 2) == 0:
        signAR = -1

    # print mu, widthL, widthR, alphaL, alphaR, nL, nR
    AL = signAL * pow(abs(nL)/abs(alphaL), nL) * np.exp(-abs(alphaL*alphaL)/2.0)
    AR = signAR * pow(abs(nR)/abs(alphaR), nR) * np.exp(-abs(alphaR*alphaR)/2.0)
    BL = nL/abs(alphaL) - abs(alphaL)
    BR = nR/abs(alphaR) - abs(alphaR)

    diff = (x[0]-mu)/width

    if diff < -alphaL:
        result = AL * pow(abs(BL-diff), -nL)
    elif -alphaL < diff < alphaR:
        result = np.exp(-0.5 * diff*diff)
    elif diff > alphaR:
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
    root_hist.GetYaxis().SetNdivisions(505)  # Unterteilung der y-Achse
    root_hist.GetYaxis().SetMaxDigits(3)  # 3 ist die einzig sinnvolle Einstellung, weil als Exponent der Zehnerpotenz nur Vielfache von 3 verwendet werden.
    root_hist.GetYaxis().SetTitleOffset(1.5)

    root_hist.Draw('HIST')
    root_hist.GetFunction(fit_function).Draw('SAME')
    legend.Draw()
    c.Print(plot_directory+filename_graphic)


if __name__ == '__main__':
    subfolder = '/generation_32'
    filename = 'histogram_files/correlator_hist_trip_32.root'
    pt_jet_range = (450, 500)
    peak_window = (0.7, 1.7)

    mtop_bw_names = ['{:.2f}'.format(elem) if elem is not None else '172.50' for elem in rew_samples]
    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    for mtop_bw, mtop_bw_name in zip(rew_samples, mtop_bw_names):
        hist = prepare_histogram(filename_root_hist=filename,
                                 hist_name='/Top-Quark/Gen-Level/weighted/correlator_hist_Gen_{:}_{:}_{:}'.format(mtop_bw, pt_jet_range[0], pt_jet_range[1]),
                                 bin_number=45)

        # Find the maximum bin
        peak_max_bin, peak_max_value = find_peak_argmax(root_hist=hist, window=peak_window)
        print(peak_max_bin, peak_max_value)

        # Perform a Gauss fit to get a rough idea of fit parameters
        gauss_fit_func, gauss_norm, gauss_mean, gauss_sigma = perform_gauss_fit(root_hist=hist, window=peak_window)
        print(gauss_norm, gauss_mean, gauss_sigma)
        draw_histogram(root_hist=hist, filename_graphic=subfolder+'/peak_fitting/correlator_gauss_fit_Gen_{:}_{:}-{:}.pdf'.format(mtop_bw_name, pt_jet_range[0], pt_jet_range[1]),
                       sample_name=mtop_bw_name, fit_label='Gauss-Fit', fit_function='gaussFitFunction')

        # Now fit with double-sided crystal ball function
        perform_crystal_ball_fit(root_hist=hist, window=peak_window, norm_gauss=gauss_norm, mean_gauss=gauss_mean, sigma_gauss=gauss_sigma)
        draw_histogram(root_hist=hist, filename_graphic=subfolder+'/peak_fitting/correlator_CB_fit_Gen_{:}_{:}-{:}.pdf'.format(mtop_bw_name, pt_jet_range[0], pt_jet_range[1]),
                       sample_name=mtop_bw_name, fit_label='Crystal-Ball-Fit', fit_function='CBFunction')

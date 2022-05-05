import ROOT
import Analysis.Tools.syncer                    # Starts syncing by itself, does not need to be called in script
from MTopCorrelations.Tools.user import plot_directory


def style_corr_hist(filename_root, hist_name, filename_graphic, ylim=(0, 0.0005), fit=False, verb=True):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hist = f.Get(hist_name)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    if fit:
        g1 = ROOT.TF1('gauss_fit', 'gaus', 1, 2.5)
        fit_result = hist.Fit(g1, 'QSR')
        fit_parameter = fit_result.Parameter(1)
        if 1 < fit_parameter < 2.5:
            print('{:.3f}'.format(fit_parameter))

    # hist.Rebin(4)
    hist.SetLineColor(ROOT.kRed)
    hist.SetTitle('')
    hist.SetLineWidth(2)
    hist.SetFillColor(ROOT.kRed)
    hist.SetLineStyle(1)
    hist.GetXaxis().SetRangeUser(0, 3)    # x axis range (also works for y axis)
    hist.GetXaxis().SetTitle("3#zeta")
    hist.GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.GetYaxis().SetRangeUser(ylim[0], ylim[1])
    hist.GetYaxis().SetTitle("Energy-weighted Triplets")
    hist.GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.Draw('HIST')

    if fit:
        g1.Draw('SAME')
        g1.SetLineColor(ROOT.kBlue)

    c.Print(plot_directory+filename_graphic)

    if verb:
        hist.GetXaxis().SetRangeUser(1, 3)
        peak_mean = hist.GetMean()
        print('{:.3f}'.format(peak_mean))


def style_top_pt_hist(filename_root, hist_name, filename_graphic):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hist = f.Get(hist_name)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    hist.SetLineColor(ROOT.kRed)
    hist.SetTitle('')
    hist.SetLineWidth(2)
    hist.SetFillColor(ROOT.kRed)
    hist.SetLineStyle(1)
    hist.GetXaxis().SetRangeUser(0, 3)    # x axis range (also works for y axis)
    hist.GetXaxis().SetTitle("p_{T,top}")
    hist.GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.GetYaxis().SetRangeUser(0, 5000000)
    hist.GetYaxis().SetTitle("Occurrence")
    hist.GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.Draw('HIST')

    c.Print(plot_directory+filename_graphic)


def style_corr_hist_unweighted(filename_root, hist_name, filename_graphic):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hist = f.Get(hist_name)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    # hist.Rebin(4)
    hist.SetLineColor(ROOT.kRed)
    hist.SetTitle('')
    hist.SetLineWidth(2)
    hist.SetFillColor(ROOT.kRed)
    hist.SetLineStyle(1)
    hist.GetXaxis().SetRangeUser(0, 3)    # x axis range (also works for y axis)
    hist.GetXaxis().SetTitle("3#zeta")
    hist.GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    # hist.GetYaxis().SetRangeUser(0, 0.002)
    hist.GetYaxis().SetTitle("Number of Triplets")
    hist.GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.Draw('HIST')

    c.Print(plot_directory+filename_graphic)


def style_numb_triplets_hist(filename_root, hist_name, filename_graphic):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hist = f.Get(hist_name)
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    # hist.Rebin(4)
    hist.SetLineColor(ROOT.kRed)
    hist.SetTitle('')
    hist.SetLineWidth(2)
    hist.SetFillColor(ROOT.kRed)
    hist.SetLineStyle(1)
    hist.GetXaxis().SetRangeUser(0, 100000)    # x axis range (also works for y axis)
    hist.GetXaxis().SetTitle(hist_name)
    hist.GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    # hist.GetYaxis().SetRangeUser(0, 0.002)
    hist.GetYaxis().SetTitle("Number of Events")
    hist.GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.Draw('HIST')

    c.Print(plot_directory+filename_graphic)


if __name__ == '__main__':
    delta_delta = 0.09

    for pt_range in ['', '_400_450', '_450_500', '_500_550', '_550_600', '_600_650', '_650_700']:
        style_corr_hist(filename_root='correlator-joined_delta_{:.2f}.root'.format(delta_delta),
                        hist_name='correlator_hist'+pt_range,
                        filename_graphic='/correlator_hist_delta_{:.2f}{:}.png'.format(delta_delta, pt_range),
                        ylim=(0, 0.0015), fit=False)
        style_top_pt_hist(filename_root='correlator-joined_delta_{:.2f}.root'.format(delta_delta),
                          hist_name='top_pt_hist'+pt_range,
                          filename_graphic='/top_pt_hist_delta_{:.2f}{:}.png'.format(delta_delta, pt_range))

    style_corr_hist_unweighted(filename_root='correlator-joined_delta_{:.2f}.root'.format(delta_delta),
                               hist_name='correlator_hist_unweighted',
                               filename_graphic='/correlator_hist_unweighted_delta_{:.2f}.png'.format(delta_delta))

    for plot_name in ['number_of_all_triplets', 'number_of_equidistant_triplets']:
        style_numb_triplets_hist(filename_root='correlator-joined_delta_{:.2f}.root'.format(delta_delta),
                                 hist_name=plot_name,
                                 filename_graphic='/'+plot_name+'_delta_{:.2f}.png'.format(delta_delta))

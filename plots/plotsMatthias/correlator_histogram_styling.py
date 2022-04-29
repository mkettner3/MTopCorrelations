import ROOT
import Analysis.Tools.syncer                    # Starts syncing by itself, does not need to be called in script
from MTopCorrelations.Tools.user import plot_directory


def style_corr_hist(filename_root, hist_name, filename_graphic, ylim=(0, 0.002)):
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
    hist.GetXaxis().SetRangeUser(0.3, 6)    # x axis range (also works for y axis)
    hist.GetXaxis().SetTitle("3#zeta")
    hist.GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hist.GetYaxis().SetRangeUser(ylim[0], ylim[1])
    hist.GetYaxis().SetTitle("Energy-weighted Triplets")
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
    hist.GetXaxis().SetRangeUser(0.3, 6)    # x axis range (also works for y axis)
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
    for pt_range in ['', '_400_450', '_450_500', '_500_550', '_550_600', '_600_650', '_650_700']:
        style_corr_hist(filename_root='correlator-joined.root',
                        hist_name='correlator_hist'+pt_range,
                        filename_graphic='/correlator_hist'+pt_range+'.png', ylim=(0, 0.00075))

    style_corr_hist_unweighted(filename_root='correlator-joined.root',
                               hist_name='correlator_hist_unweighted',
                               filename_graphic='/correlator_hist_unweighted.png')

    for plot_name in ['number_of_all_triplets', 'number_of_equidistant_triplets']:
        style_numb_triplets_hist(filename_root='correlator-joined.root',
                                 hist_name=plot_name,
                                 filename_graphic='/'+plot_name+'.png')
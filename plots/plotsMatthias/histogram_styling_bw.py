import ROOT
from MTopCorrelations.Tools.user import plot_directory
from calc_triplet_and_hist import pt_jet_ranges


def style_corr_hist(filename_root, hist_name, sample_names, filename_graphic, ylim=(0, 0.0005), verb=True):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hists = [f.Get(hist_name.replace('$', sample_name)) for sample_name in sample_names]
    for i in range(len(hists)):
        hists[i].SetDirectory(ROOT.nullptr)
    f.Close()
    c = ROOT.TCanvas('c', 'c', 600, 600)
    legend = ROOT.TLegend(0.75, 0.6, 0.94, 0.89)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    if len(hists) == 3:
        line_colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kRed]
    elif len(hists) == 5:
        line_colors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kYellow]
    elif len(hists) == 9:
        line_colors = list(range(1, 10))
    else:
        raise RuntimeError('Please specify the line colors in style_corr_hist()!')

    for hist, line_color, sample_name in zip(hists, line_colors, sample_names):
        if sample_name == 'None':
            sample_name = '172.5'
        hist.Rebin(20)
        hist.SetLineColor(line_color)
        hist.SetTitle('')
        hist.SetLineWidth(2)
        hist.SetLineStyle(1)
        legend.AddEntry(hist, sample_name[-5:], 'l')
        if hist_name[-6:] != 'abscou':
            hist.Scale(1/hist.Integral())
    hists[0].GetXaxis().SetRangeUser(0, 3)    # x-axis range (also works for y-axis)
    hists[0].GetXaxis().SetTitle("3#zeta")
    hists[0].GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hists[0].GetYaxis().SetRangeUser(ylim[0], ylim[1])
    hists[0].GetYaxis().SetTitle("Energy-weighted Triplets")
    hists[0].GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse

    hists[0].Draw('HIST')
    for i in range(1, len(hists)):
        hists[i].Draw('HIST SAME')

    legend.Draw()
    c.Print(plot_directory+filename_graphic)

    if verb:
        hists[0].GetXaxis().SetRangeUser(1, 3)
        peak_mean = hists[0].GetMean()
        print('{:.3f}'.format(peak_mean))


def style_varied_hist(filename_root, hist_name, varied_hist_name, var_factors, filename_graphic, ylim=(0, 0.0005), verb=True):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hists = {}
    hists[int(len(var_factors)/2)] = f.Get(hist_name)         # Only works if the length of var_factors is even
    for i in range(int(len(var_factors)/2)):
        hists[i] = f.Get(varied_hist_name.replace('$', var_factors[i]))
    for i in range(int(len(var_factors)/2+1), len(var_factors)+1):
        hists[i] = f.Get(varied_hist_name.replace('$', var_factors[i-1]))

    for i in range(len(hists)):
        hists[i].SetDirectory(ROOT.nullptr)
    f.Close()
    c = ROOT.TCanvas('c', 'c', 600, 600)
    legend = ROOT.TLegend(0.75, 0.7, 0.94, 0.89)
    for i, label in enumerate([var_factors[j] for j in range(int(len(var_factors)/2))] + ['original'] + [var_factors[j] for j in range(int(len(var_factors)/2), len(var_factors))]):
        legend.AddEntry(hists[i], label, 'l')
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    if len(hists) == 3:
        line_colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kRed]
    elif len(hists) == 5:
        line_colors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kYellow]
    elif len(hists) == 7:
        line_colors = [ROOT.kMagenta, ROOT.kViolet, ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kOrange, ROOT.kYellow]
    elif len(hists) == 9:
        line_colors = list(range(1, 10))
    else:
        raise RuntimeError('Please specify the line colors in style_corr_hist()!')

    for hist, line_color in zip(hists.values(), line_colors):
        hist.Rebin(20)
        hist.SetLineColor(line_color)
        hist.SetTitle('')
        hist.SetLineWidth(2)
        hist.SetLineStyle(1)
        if hist_name[-6:] != 'abscou':
            hist.Scale(1/hist.Integral())
    hists[0].GetXaxis().SetRangeUser(0, 3)    # x-axis range (also works for y-axis)
    hists[0].GetXaxis().SetTitle("3#zeta")
    hists[0].GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hists[0].GetYaxis().SetRangeUser(ylim[0], ylim[1])
    hists[0].GetYaxis().SetTitle("Energy-weighted Triplets")
    hists[0].GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse

    hists[0].Draw('HIST')
    for i in range(1, len(hists)):
        hists[i].Draw('HIST SAME')

    legend.Draw()
    c.Print(plot_directory+filename_graphic)

    if verb:
        hists[0].GetXaxis().SetRangeUser(1, 3)
        peak_mean = hists[0].GetMean()
        print('{:.3f}'.format(peak_mean))


def style_jet_hist(filename_root, sample_names, hist_name, filename_graphic, xlim=(380, 730), ylim=(0, 0.0005), title='Hadronic Top-Jet-Mass', verb=True):
    ROOT.gStyle.SetLegendBorderSize(0)  # No border for legend
    ROOT.gStyle.SetPadTickX(1)          # Axis ticks on top
    ROOT.gStyle.SetPadTickY(1)          # Axis ticks right
    ROOT.gStyle.SetOptStat(0)           # Do not display stat box

    f = ROOT.TFile(filename_root)
    hists = [f.Get(hist_name.replace('$', sample_name)) for sample_name in sample_names]
    for i in range(len(hists)):
        hists[i].SetDirectory(ROOT.nullptr)
    f.Close()
    c = ROOT.TCanvas('c', 'c', 600, 600)
    ROOT.gPad.SetLeftMargin(0.19)
    ROOT.gPad.SetBottomMargin(0.2)

    if len(hists) == 3:
        line_colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kRed]
    elif len(hists) == 2:
        line_colors = [ROOT.kGreen, ROOT.kRed]
    elif len(hists) == 4:
        line_colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kYellow]
    elif len(hists) == 5:
        line_colors = [ROOT.kMagenta, ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kYellow]
    elif len(hists) == 7:
        line_colors = [ROOT.kMagenta, ROOT.kViolet, ROOT.kBlue, ROOT.kGreen, ROOT.kRed, ROOT.kOrange, ROOT.kYellow]
    elif len(hists) == 9:
        line_colors = list(range(1, 10))
    else:
        raise RuntimeError('Please specify the line colors in style_corr_hist()!')

    for hist, line_color, sample_name in zip(hists, line_colors, sample_names):
        hist.Rebin(2)
        hist.SetLineColor(line_color)
        hist.SetTitle('')
        hist.SetLineWidth(2)
        hist.SetLineStyle(1)
    hists[0].GetXaxis().SetRangeUser(xlim[0], xlim[1])
    hists[0].GetXaxis().SetTitle(title)
    hists[0].GetXaxis().SetNdivisions(505)      # Unterteilung der x-Achse
    hists[0].GetYaxis().SetRangeUser(ylim[0], ylim[1])
    hists[0].GetYaxis().SetTitle("Number of Events")
    hists[0].GetYaxis().SetNdivisions(505)      # Unterteilung der x-Achse

    hists[0].Draw('HIST')
    for i in range(1, len(hists)):
        hists[i].Draw('HIST SAME')

    c.Print(plot_directory+filename_graphic)

    if verb:
        hists[0].GetXaxis().SetRangeUser(1, 3)
        peak_mean = hists[0].GetMean()
        print('{:.3f}'.format(peak_mean))


if __name__ == '__main__':
    subfolder = '/generation_31'
    filename = 'histogram_files/correlator_hist_trip_31.root'
    sample_names = ['171.5', '171.75', '172.0', '172.25', 'None', '172.75', '173.0', '173.25', '173.5']

    ROOT.gROOT.SetBatch(ROOT.kTRUE)             # Prevent graphical display for every c.Print() statement

    for level in ['Gen', 'PF']:
        for pt_range in pt_jet_ranges:
            style_corr_hist(filename_root=filename,
                            hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_$_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
                            sample_names=sample_names,
                            filename_graphic=subfolder+'/correlator_hist_{:}_{:}-{:}.png'.format(level, pt_range[0], pt_range[1]),
                            ylim=(0, 0.006), verb=False)

            style_corr_hist(filename_root=filename,
                            hist_name='/Top-Quark/'+level+'-Level/absolute/correlator_hist_{:}_$_{:}_{:}_abscou'.format(level, pt_range[0], pt_range[1]),
                            sample_names=sample_names,
                            filename_graphic=subfolder+'/correlator_hist_{:}_{:}-{:}_abs.png'.format(level, pt_range[0], pt_range[1]),
                            ylim=(0, 500000), verb=False)

            style_corr_hist(filename_root=filename,
                            hist_name='/W-Boson/'+level+'-Level/weighted/correlator_hist_W_{:}_$_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
                            sample_names=sample_names,
                            filename_graphic=subfolder+'/correlator_hist_W_{:}_{:}-{:}.png'.format(level, pt_range[0], pt_range[1]),
                            ylim=(0, 0.08), verb=False)

            style_corr_hist(filename_root=filename,
                            hist_name='/W-Boson/'+level+'-Level/absolute/correlator_hist_W_{:}_$_{:}_{:}_abscou'.format(level, pt_range[0], pt_range[1]),
                            sample_names=sample_names,
                            filename_graphic=subfolder+'/correlator_hist_W_{:}_{:}-{:}_abs.png'.format(level, pt_range[0], pt_range[1]),
                            ylim=(0, 400000), verb=False)

            style_varied_hist(filename_root=filename,
                              hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_None_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
                              varied_hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_jet_$_{:}_None_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
                              var_factors=['{:.2f}'.format(value) for value in [1.1, 1.02, 0.98, 0.9]],
                              filename_graphic=subfolder+'/correlator_hist_varied_jet_{:}_{:}-{:}.png'.format(level, pt_range[0], pt_range[1]),
                              ylim=(0, 0.006), verb=False)

            style_varied_hist(filename_root=filename,
                              hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_None_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
                              varied_hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_cons_pt_$_{:}_None_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
                              var_factors=['{:.2f}'.format(value) for value in [-2, -1, -0.5, 0.5, 1, 2]],
                              filename_graphic=subfolder+'/correlator_hist_varied_cons_pt_{:}_{:}-{:}.png'.format(level, pt_range[0], pt_range[1]),
                              ylim=(0, 0.006), verb=False)

            # style_varied_hist(filename_root=filename,
            #                   hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_{:}_None_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
            #                   varied_hist_name='/Top-Quark/'+level+'-Level/weighted/correlator_hist_varied_cons_eta_phi_$_{:}_None_{:}_{:}'.format(level, pt_range[0], pt_range[1]),
            #                   var_factors=['1.02', '0.98'],
            #                   filename_graphic=subfolder+'/correlator_hist_varied_cons_eta_phi_{:}_{:}-{:}.png'.format(level, pt_range[0], pt_range[1]),
            #                   ylim=(0, 0.02), verb=False)

        style_jet_hist(filename_root=filename,
                       hist_name='/Others/'+level+'-Level/hadronic_top_jet_pt_hist_{:}_$'.format(level),
                       sample_names=sample_names,
                       filename_graphic=subfolder+'/hadronic_top_jet_pt_hist_{:}.png'.format(level),
                       xlim=(380, 730), ylim=(0, 200), title='Hadronic Top-Jet-p_{T}', verb=False)

        style_jet_hist(filename_root=filename,
                       hist_name='/Others/'+level+'-Level/hadronic_top_jet_mass_hist_{:}_$'.format(level),
                       sample_names=sample_names,
                       filename_graphic=subfolder+'/hadronic_top_jet_mass_hist_{:}.png'.format(level),
                       xlim=(160, 195), ylim=(0, 800), title='Hadronic Top-Jet-Mass', verb=False)

        style_jet_hist(filename_root=filename,
                       hist_name='/Others/'+level+'-Level/hadronic_top_jet_mass_hist_{:}_MC_$'.format(level),
                       sample_names=['TTbar_1', 'TTbar_2', 'TTbar_4', 'TTbar_5'],
                       filename_graphic=subfolder+'/hadronic_top_jet_mass_hist_{:}_MC.png'.format(level),
                       xlim=(160, 195), ylim=(0, 800), title='Hadronic Top-Jet-Mass (Monte-Carlo)', verb=False)

        style_jet_hist(filename_root=filename,
                       hist_name='/Others/'+level+'-Level/hadronic_top_jet_mass_hist_{:}_$'.format(level),
                       sample_names=['173.5', 'MC_TTbar_4'],
                       filename_graphic=subfolder+'/hadronic_top_jet_mass_hist_{:}_MC_comparison.png'.format(level),
                       xlim=(160, 195), ylim=(0, 800), title='Hadronic Top-Jet-Mass (Comparison)', verb=False)

        style_jet_hist(filename_root=filename,
                       hist_name='/Others/'+level+'-Level/hadronic_top_mass_hist_{:}_$'.format(level),
                       sample_names=sample_names,
                       filename_graphic=subfolder+'/hadronic_top_mass_hist_{:}.png'.format(level),
                       xlim=(170, 175), ylim=(0, 600), title='Hadronic Top-Mass', verb=False)

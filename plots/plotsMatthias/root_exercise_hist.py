# -*- coding: utf-8 -*-

import numpy as np
from ROOT import TH1F, TCanvas, TFile, kRed

hist = TH1F('exercise_hist', 'Exercise Histogram', 10, 0, 100)
hist.Fill(11.3)
hist.Fill(25.4)
hist.Fill(18.1)
for number in np.square(np.arange(10)):
    hist.Fill(number)
print('Histrogram was created and filled.')

hist_mean = hist.GetMean()
hist_rms = hist.GetRMS()

print('Mean value of the histogram: {:.2f}'.format(hist_mean))
print('RMS value of the histogram: {:.2f}'.format(hist_rms))

hist_integral = hist.Integral()
max_bin = hist.GetMaximumBin()
max_bin_content = hist.GetBinContent(max_bin)

hist.GetYaxis().SetTitle('entries')
hist.SetLineColor(kRed)

c = TCanvas("canvas", "canvas", 800, 600)
hist.Draw()
c.Print('exercise_hist.png')

f = TFile('tests.root', 'RECREATE')
f.cd()
hist.Write('funny_histogram')
f.Close()

# -*- coding: utf-8 -*-

"""
Compute the statistical uncertainty as expected from data
"""

from calc_chi2 import calc_norm_cov_matrix


pt_jet_lowest = 400
pt_jet_highest = 700
pt_jet_step = 50
pt_jet_ranges = zip(range(pt_jet_lowest, pt_jet_highest, pt_jet_step),
                    range(pt_jet_lowest+50, pt_jet_highest+50, pt_jet_step))

subfolder = '/uncertainty_from_data'
filename = 'histogram_files/correlator_hist_trip_8.root'
sample_names = ['TTbar_169p5', 'TTbar_171p5', 'TTbar_172p5', 'TTbar_173p5', 'TTbar_175p5']

matrices_norm = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
root_hist = [[[None for _ in range(len(pt_jet_ranges))] for _ in range(len(sample_names))] for _ in range(2)]
hist_axis_range = (0, 3)

for g, level in enumerate(['Gen', 'PF']):
    for h, sample_name in enumerate(sample_names):
        print('Working on sample "{}" ...'.format(sample_name))
        for k, pt_range in enumerate(pt_jet_ranges):
            (matrices_norm[g][h][k],
             root_hist[g][h][k],
             hist_axis_range) = calc_norm_cov_matrix(filename_root_hist=filename,
                                                     hist_name='/Top-Quark/'+level+'-Level/absolute/correlator_hist_{:}_{:}_{:}_{:}_abscou'.format(level, sample_name,
                                                                                                        pt_range[0], pt_range[1]),
                                                     plot_matrix=True,
                                                     id_level=level, id_sample=sample_name, id_range=pt_range,
                                                     absolute_hist=True)

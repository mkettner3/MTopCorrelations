# -*- coding: utf-8 -*-

"""
Functions to enable conversions from ROOT datatypes to numpy datatypes
"""

from typing import Any
import numpy as np


def root_hist_to_numpy(root_hist):
    # type: (Any) -> np.ndarray

    if root_hist.Class_Name()[3] == 'F':
        numpy_type = np.single
    elif root_hist.Class_Name()[3] == 'D':
        numpy_type = np.double
    elif root_hist.Class_Name()[3] == 'I':
        numpy_type = np.intc
    elif root_hist.Class_Name()[3] == 'S':
        numpy_type = np.short
    else:
        raise ValueError('Only the root datatypes F, D, I and S are supported by this function!')

    if root_hist.Class_Name()[:3] == 'TH1':
        bins = root_hist.GetNbinsX()
        array = np.zeros(bins, dtype=numpy_type)

        for numb_bin in range(bins):
            array[numb_bin] = root_hist.GetBinContent(numb_bin+1)

    elif root_hist.Class_Name()[:3] == 'TH2':
        x_bins = root_hist.GetNbinsX()
        y_bins = root_hist.GetNbinsY()
        array = np.zeros((x_bins, y_bins), dtype=numpy_type)

        for y_bin in range(y_bins):
            for x_bin in range(x_bins):
                array[x_bin, y_bin] = root_hist.GetBinContent(x_bin+1, y_bin+1)

    else:
        raise ValueError('The object "root_hist" must be a 1 or 2-dimensional histogram!')

    return array

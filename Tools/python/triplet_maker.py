# -*- coding: utf-8 -*-

"""
Function to produce triplets.
"""

import numpy as np


def make_triplets(jet_pt, particle_vectors, n=2, top_data=True, w_data=True, pt_value=True):
    # type: (float, list, int, bool, bool, bool) -> np.ndarray
    """
    Method to produce triplets from the particle vectors.

    :param jet_pt: float; The Pt of the jet from the top quark.
    :param particle_vectors: list of ROOT.TLorentzVector(); The lorentz-vectors of the particles inside the
        mentioned jet of one event.
    :param n: integer; Order of energy weighting. Default is quadratic order (n=2).
    :param top_data: boolean; Decides if the maximum of delta_zeta is returned within the numpy array.
    :param w_data: boolean; Decides if the delta_zeta between legs and shortest side are returned within the
        numpy array.
    :param pt_value: boolean; Decides if the jet_pt value is returned within the numpy array.

    :return: two-dimensional numpy array; Rows are triplets and Columns are [3*Zeta, Weight (w), maximum of
        delta_zeta, the delta_zeta between the legs of an isosceles triangle, the length of the shortest triangle
        side, the pt of the surrounding jet]
    """
    w = []
    zeta_value = []

    for i in range(len(particle_vectors)):
        for j in range(i+1, len(particle_vectors)):
            for k in range(j+1, len(particle_vectors)):
                zeta_value.append([particle_vectors[i].DeltaR(particle_vectors[j]),
                                   particle_vectors[i].DeltaR(particle_vectors[k]),
                                   particle_vectors[j].DeltaR(particle_vectors[k])])

                w.append((particle_vectors[i].Pt() * particle_vectors[j].Pt() * particle_vectors[k].Pt())**n / ((jet_pt**3)**n))

    zeta_value = np.array(zeta_value, dtype=np.float32)
    zeta_value.sort(axis=1)
    three_zeta = zeta_value.sum(axis=1)

    triplet_param_list = [three_zeta, np.array(w, dtype=np.float32)]

    if top_data:
        triplet_param_list.append(zeta_value[:, 0] - zeta_value[:, 2])          # max_delta_zeta

    if w_data:
        triplet_param_list.append(zeta_value[:, 0] - zeta_value[:, 1])          # delta_legs
        triplet_param_list.append(zeta_value[:, 2])                             # shortest_side

    if pt_value:
        triplet_param_list.append(np.full((len(w)), jet_pt, dtype=np.float32))

    triplet_params = np.stack(triplet_param_list, axis=1)

    return triplet_params

# -*- coding: utf-8 -*-

"""
Function to produce triplets.
"""

import numpy as np


class Triplets(object):

    @staticmethod
    def make(jet_pt, particle_vectors, n=2):
        # type: (float, list, int) -> np.ndarray
        """
        Method to produce triplets from the particle vectors.

        :param jet_pt: float; The Pt of the jet from the top quark.
        :param particle_vectors: list of ROOT.TLorentzVector(); The lorentz-vectors of the particles inside the
            mentioned jet of one event.
        :param n: integer; Order of energy weighting. Default is quadratic order (n=2).

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

        max_delta_zeta = zeta_value[:, 0] - zeta_value[:, 2]
        delta_legs = zeta_value[:, 0] - zeta_value[:, 1]
        shortest_side = zeta_value[:, 2]

        triplet_params = np.stack([three_zeta, np.array(w, dtype=np.float32), max_delta_zeta,
                                   delta_legs, shortest_side, np.full((len(w)), jet_pt, dtype=np.float32)], axis=1)

        return triplet_params

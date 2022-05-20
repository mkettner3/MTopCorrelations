# -*- coding: utf-8 -*-

"""
Function to produce triplets.
"""

import numpy as np


class Triplets(object):

    @staticmethod
    def make(jet_pt: float, particle_vectors: list) -> np.ndarray:
        """
        Method to produce triplets from the particle vectors.

        :param jet_pt: float; The Pt of the jet from the top quark.
        :param particle_vectors: list of ROOT.TLorentzVector(); The lorentz-vectors of the particles inside the
            mentioned jet of one event.

        :return: two-dimensional numpy array; Row are triplets and columns are [Zeta, Weight (w), maximum of delta_zeta,
            the delta_zeta between the legs of an isosceles triangle, the length of the shortest triangle side]
        """
        triplet_params = []

        for i in range(len(particle_vectors)):
            for j in range(i+1, len(particle_vectors)):
                for k in range(j+1, len(particle_vectors)):
                    zeta_0_1 = particle_vectors[i].DeltaR(particle_vectors[j])
                    zeta_0_2 = particle_vectors[i].DeltaR(particle_vectors[k])
                    zeta_1_2 = particle_vectors[j].DeltaR(particle_vectors[k])

                    w = (particle_vectors[i].Pt() * particle_vectors[j].Pt() * particle_vectors[k].Pt())**2 / (jet_pt**6)
                                                                                # 6 = 2*3
                    zeta = (zeta_0_1 + zeta_0_2 + zeta_1_2) / 3

                    max_delta_zeta = max(abs(zeta_0_1-zeta_0_2), abs(zeta_1_2-zeta_0_1), abs(zeta_0_2-zeta_1_2))
                    zeta_list = sorted([zeta_0_1, zeta_0_2, zeta_1_2], reverse=True)
                    delta_legs = abs(zeta_list[0] - zeta_list[1])
                    shortest_side = min(zeta_0_1, zeta_0_2, zeta_1_2)

                    triplet_params.append([zeta, w, max_delta_zeta, delta_legs, shortest_side])

        return np.array(triplet_params, dtype=np.float32)

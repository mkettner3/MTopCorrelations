# -*- coding: utf-8 -*-

"""
Obtain the constituents of a jet.
"""

from ROOT import TLorentzVector
from typing import Any


class JetConstituents(object):

    @staticmethod
    def get(event, index, max_numb_of_cons=None):
        # type: (Any, int, int) -> list
        """
        Method to obtain the lorentz-vectors of the constituents of a jet.

        :param event: A ROOT event.
        :param index: int; The index of the respective jet.
        :param max_numb_of_cons: int; The maximum number of jet constituents, which should be returned.

        :return: list; A list of lorentz-vectors of the constituents of the jet.
        """
        constituents = []
        for i in range(event.nGenJetAK8_cons):
            if event.GenJetAK8_cons_jetIndex[i] == index:
                particle = TLorentzVector()
                particle.SetPtEtaPhiM(event.GenJetAK8_cons_pt[i], event.GenJetAK8_cons_eta[i],
                                      event.GenJetAK8_cons_phi[i], event.GenJetAK8_cons_mass[i])
                constituents.append(particle)

        def get_pt(part):
            return part.Pt()

        constituents.sort(key=get_pt, reverse=True)
        if isinstance(max_numb_of_cons, int):
            if len(constituents) > max_numb_of_cons:
                constituents = constituents[:max_numb_of_cons]

        return constituents

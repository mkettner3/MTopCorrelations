# -*- coding: utf-8 -*-

"""
Obtain the constituents of a jet.
"""

from ROOT import TLorentzVector


class JetConstituents(object):

    @staticmethod
    def get(event, index) -> list:
        """
        Method to obtain the lorentz-vectors of the constituents of a jet.

        :param event: float; A ROOT event.
        :param index: int; The index of the respective jet.

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
        return constituents

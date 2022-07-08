# -*- coding: utf-8 -*-

"""
Obtain the constituents of a jet.
"""

from ROOT import TLorentzVector
from typing import Any


def get_jet_constituents(event, level, index, max_numb_of_cons=None):
    # type: (Any, str, int, int) -> list
    """
    Method to obtain the lorentz-vectors of the constituents of a jet.

    :param event: A ROOT event.
    :param level: str; 'Gen' or 'PF' determines if the particle level or the detector level is used.
    :param index: int; The index of the respective jet.
    :param max_numb_of_cons: int; The maximum number of jet constituents, which should be returned.

    :return: list; A list of lorentz-vectors of the constituents of the jet.
    """
    if level == 'Gen':
        event_para = {'Jet_cons': event.nGenJetAK8_cons, 'Jet_cons_jetIndex': event.GenJetAK8_cons_jetIndex,
                      'Jet_cons_pt': event.GenJetAK8_cons_pt, 'Jet_cons_eta': event.GenJetAK8_cons_eta,
                      'Jet_cons_phi': event.GenJetAK8_cons_phi, 'Jet_cons_mass': event.GenJetAK8_cons_mass}
    elif level == 'PF':
        event_para = {'Jet_cons': event.nPFJetAK8_cons, 'Jet_cons_jetIndex': event.PFJetAK8_cons_jetIndex,
                      'Jet_cons_pt': event.PFJetAK8_cons_pt, 'Jet_cons_eta': event.PFJetAK8_cons_eta,
                      'Jet_cons_phi': event.PFJetAK8_cons_phi, 'Jet_cons_mass': event.PFJetAK8_cons_mass}
    else:
        raise ValueError('The parameter "level" must be either "Gen" or "PF"!')

    constituents = []
    for i in range(event_para['Jet_cons']):
        if event_para['Jet_cons_jetIndex'][i] == index:
            particle = TLorentzVector()
            particle.SetPtEtaPhiM(event_para['Jet_cons_pt'][i], event_para['Jet_cons_eta'][i],
                                  event_para['Jet_cons_phi'][i], event_para['Jet_cons_mass'][i])
            constituents.append(particle)

    def get_pt(part):
        return part.Pt()

    constituents.sort(key=get_pt, reverse=True)
    if isinstance(max_numb_of_cons, int):
        if len(constituents) > max_numb_of_cons:
            constituents = constituents[:max_numb_of_cons]

    return constituents

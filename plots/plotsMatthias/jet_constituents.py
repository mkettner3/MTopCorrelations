# -*- coding: utf-8 -*-

"""
Obtain the constituents of a jet.
"""

from ROOT import TLorentzVector
from typing import Any


def get_jet_constituents(event, level, index, particle_types='all', max_numb_of_cons=None):
    # type: (Any, str, int, str, int) -> list
    """
    Method to obtain the lorentz-vectors of the constituents of a jet.

    :param event: A ROOT event.
    :param level: str; 'Gen' or 'PF' determines if the particle level or the detector level is used.
    :param index: int; The index of the respective jet.
    :param particle_types: string; Determines which particle types are returned. It can be 'all', 'charged_par' or
        'charged_had', to return either all constituents, all charged particles or only charged hadrons.
    :param max_numb_of_cons: int; The maximum number of jet constituents, which should be returned.

    :return: list; A list of lorentz-vectors of the constituents of the jet.
    """
    if level == 'Gen':
        event_para = {'Jet_cons': event.nGenJetAK8_cons, 'Jet_cons_jetIndex': event.GenJetAK8_cons_jetIndex,
                      'Jet_cons_pt': event.GenJetAK8_cons_pt, 'Jet_cons_eta': event.GenJetAK8_cons_eta,
                      'Jet_cons_phi': event.GenJetAK8_cons_phi, 'Jet_cons_mass': event.GenJetAK8_cons_mass,
                      'Jet_cons_pdgId': event.GenJetAK8_cons_pdgId}
    elif level == 'PF':
        event_para = {'Jet_cons': event.nPFJetAK8_cons, 'Jet_cons_jetIndex': event.PFJetAK8_cons_jetIndex,
                      'Jet_cons_pt': event.PFJetAK8_cons_pt, 'Jet_cons_eta': event.PFJetAK8_cons_eta,
                      'Jet_cons_phi': event.PFJetAK8_cons_phi, 'Jet_cons_mass': event.PFJetAK8_cons_mass,
                      'Jet_cons_pdgId': event.PFJetAK8_cons_pdgId}
    else:
        raise ValueError('The parameter "level" must be either "Gen" or "PF"!')

    constituents = []
    for i in range(event_para['Jet_cons']):
        if event_para['Jet_cons_jetIndex'][i] == index:
            if particle_types == 'all' or \
                    (particle_types == 'charged_par' and abs(event_para['Jet_cons_pdgId'][i]) in [211, 13, 11, 1, 321, 2212, 3222, 3112, 3312, 3334]) or \
                    (particle_types == 'charged_had' and abs(event_para['Jet_cons_pdgId'][i]) in [211, 321, 2212, 3222, 3112, 3312, 3334]):
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

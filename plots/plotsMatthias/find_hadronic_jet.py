# -*- coding: utf-8 -*-

"""
Determine the jet from the hadronic top quark.
"""

from ROOT import TLorentzVector
import numpy as np


def find_hadronic_jet(event, level, merge_tolerance=0.8, jet_pt_min=400):
    if level == 'Gen':
        event_para = {'nJet': event.nGenJetAK8, 'Jet_pt': event.GenJetAK8_pt, 'Jet_eta': event.GenJetAK8_eta,
                      'Jet_phi': event.GenJetAK8_phi, 'Jet_mass': event.GenJetAK8_mass}
    elif level == 'PF':
        event_para = {'nJet': event.nPFJetAK8, 'Jet_pt': event.PFJetAK8_pt, 'Jet_eta': event.PFJetAK8_eta,
                      'Jet_phi': event.PFJetAK8_phi, 'Jet_mass': event.PFJetAK8_mass}
    else:
        raise ValueError('The parameter "level" must be either "Gen" or "PF"!')

    top_vec, anti_top_vec = TLorentzVector(), TLorentzVector()
    quark_vec, anti_quark_vec, bottom_vec = TLorentzVector(), TLorentzVector(), TLorentzVector()
    for i in range(event.nGenPart):
        GenPart_properties = (event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
        if event.GenPart_pdgId[i] == 6:
            top_vec.SetPtEtaPhiM(*GenPart_properties)
        elif event.GenPart_pdgId[i] == -6:
            anti_top_vec.SetPtEtaPhiM(*GenPart_properties)
        elif event.GenPart_pdgId[i] in range(1, 7) and abs(event.GenPart_grmompdgId[i]) == 6:
            quark_vec.SetPtEtaPhiM(*GenPart_properties)
        elif event.GenPart_pdgId[i] in range(-1, -7, -1) and abs(event.GenPart_grmompdgId[i]) == 6:
            anti_quark_vec.SetPtEtaPhiM(*GenPart_properties)
        elif event.GenPart_pdgId[i] == 5 and abs(event.GenPart_mompdgId[i]) == 6:
            bottom_vec.SetPtEtaPhiM(*GenPart_properties)

    top_lep, top_had = TLorentzVector(), TLorentzVector()
    for i in range(event.nGenPart):
        if abs(event.GenPart_pdgId[i]) in [11, 13, 15]:
            if event.GenPart_grmompdgId[i] == 6:
                # top_lep = top_vec
                top_had = anti_top_vec
            elif event.GenPart_grmompdgId[i] == -6:
                # top_lep = anti_top_vec
                top_had = top_vec
            break

    jets = []
    for i in range(event_para['nJet']):
        jets.append(TLorentzVector())
        jets[-1].SetPtEtaPhiM(event_para['Jet_pt'][i], event_para['Jet_eta'][i],
                              event_para['Jet_phi'][i], event_para['Jet_mass'][i])
    deltas_top = [jets[j].DeltaR(top_had) for j in range(event_para['nJet'])]
    delta_top_min = min(deltas_top)
    hadronic_jet_idx = np.argmin(deltas_top)
    hadronic_jet_pt = jets[hadronic_jet_idx].Pt()
    hadronic_jet_mass = jets[hadronic_jet_idx].M()

    delta_q = jets[hadronic_jet_idx].DeltaR(quark_vec)
    delta_aq = jets[hadronic_jet_idx].DeltaR(anti_quark_vec)
    delta_b = jets[hadronic_jet_idx].DeltaR(bottom_vec)

    if not (delta_top_min < merge_tolerance and delta_q < merge_tolerance and
            delta_aq < merge_tolerance and delta_b < merge_tolerance and hadronic_jet_pt > jet_pt_min):
        hadronic_jet_idx, hadronic_jet_pt, hadronic_jet_mass = None, None, None

    return hadronic_jet_idx, hadronic_jet_pt, hadronic_jet_mass, top_had.M()

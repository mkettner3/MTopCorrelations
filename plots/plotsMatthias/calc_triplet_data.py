# -*- coding: utf-8 -*-

"""
Script to produce numpy arrays of triplet parameters.
"""

import numpy as np
import time
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
from triplet_maker import make_triplets
from jet_constituents import get_jet_constituents
from ROOT import TLorentzVector, TNtuple, TFile
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse


# The function find_hadronic_jet() can be imported from find_hadronic_jet.py
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
        if event.GenPart_pdgId[i] == 6:
            top_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
        elif event.GenPart_pdgId[i] == -6:
            anti_top_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i],
                                      event.GenPart_m[i])
        elif event.GenPart_pdgId[i] in range(1, 7) and abs(event.GenPart_grmompdgId[i]) == 6:
            quark_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
        elif event.GenPart_pdgId[i] in range(-1, -7, -1) and abs(event.GenPart_grmompdgId[i]) == 6:
            anti_quark_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i],
                                        event.GenPart_m[i])
        elif event.GenPart_pdgId[i] == 5 and abs(event.GenPart_mompdgId[i]) == 6:
            bottom_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])

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

    return hadronic_jet_idx, hadronic_jet_pt, hadronic_jet_mass


def calc_triplet_data(sample):
    read_variables = [
        "nGenPart/I",
        "GenPart[pt/F,eta/F,phi/F,m/F,pdgId/I,mompdgId/I,grmompdgId/I]",
        "nGenJetAK8/I",
        "GenJetAK8[pt/F,eta/F,phi/F,mass/F]",
        "nGenJetAK8_cons/I",
        VectorTreeVariable.fromString("GenJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000)]

    triplets = []

    r = sample.treeReader(variables=read_variables, selectionString="GenJetAK8_pt>400")

    r.start()
    while r.run():                                                              # Event-Loop
        nearest_jet_idx, nearest_jet_pt, _ = find_hadronic_jet(r.event, level='Gen')
        jet_constituents = get_jet_constituents(event=r.event, level='Gen', index=nearest_jet_idx, max_numb_of_cons=50)

        if len(jet_constituents) > 0:
            triplets.append(make_triplets(jet_pt=nearest_jet_pt, particle_vectors=jet_constituents))

    return np.concatenate(triplets, axis=0)


def save_triplets_to_hdf5_file(triplets_data, filename, calc_time):
    # type: (np.ndarray, str, float) -> None

    import h5py

    with h5py.File(filename, 'a') as f:
        ds = f.create_dataset(name='triplet-parameters', dtype='float32', data=triplets_data)
        description = 'Rows are triplets and Columns are [3*Zeta, Weight (w), maximum of delta_zeta, the delta_zeta between the legs of an isosceles triangle, the length of the shortest triangle side, the pt of the surrounding jet]'
        ds.attrs.create(name='description', data=description)
        ds.attrs.create(name='computing time', data=calc_time, dtype='float32')


def save_triplets_to_root_file(triplets_data, filename):
    # type: (np.ndarray, str) -> None

    outfile = TFile(filename, 'RECREATE')
    outfile.cd()

    ds = TNtuple('triplet-parameters', 'triplet-parameters',
                 'three_zeta:w:max_delta_zeta:delta_legs:shortest_side:jet_pt')
    for row in triplets_data:
        ds.Fill(row)

    outfile.Close()


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    args = argParser.parse_args()

    mc_ = [UL2018.TTbar_1]            # + [UL2018.TTbar_2, UL2018.TTbar_3]
    samples = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_]

    for sample in samples:
        start = time.time()
        all_triplets = calc_triplet_data(sample=sample)
        end = time.time()
        save_triplets_to_hdf5_file(triplets_data=all_triplets,
                                   filename='triplet_files/EWC_triplets_{:}_{:02}.h5'.format(sample.name[:11], args.job),
                                   calc_time=end-start)

    # print('Executing calc_triplet_data.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))

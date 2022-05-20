# -*- coding: utf-8 -*-

"""
Script to produce numpy arrays of triplet parameters.
"""

# Define the MC samples
import numpy as np
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
from MTopCorrelations.Tools.python.triplet_maker import Triplets
from MTopCorrelations.Tools.python.jet_constituents import JetConstituents
from ROOT import TLorentzVector
import argparse


def calc_triplet_data():
    mc_ = [UL2018.TTbar_1]
    mc_dev_ = [UL2018.TTbar_2, UL2018.TTbar_3]
    mc = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_]
    mc_dev = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc_dev_]
    mc_all = mc + mc_dev

    triplets = []
    for event_loop:             # ToDo: Find out how an event loop works in ROOT
                                # ToDo: Find out how to read the necessary data efficiently from our sample files.
        top_vec, anti_top_vec = TLorentzVector(), TLorentzVector()
        for i in range(event.nGenPart):
            if event.GenPart_pdgId[i] == 6:
                top_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])
            elif event.GenPart_pdgId[i] == -6:
                anti_top_vec.SetPtEtaPhiM(event.GenPart_pt[i], event.GenPart_eta[i], event.GenPart_phi[i], event.GenPart_m[i])

        top_lep, top_had = TLorentzVector(), TLorentzVector()
        for i in range(event.nGenPart):
            if abs(event.GenPart_pdgId[i]) in [11, 13, 15]:
                if event.GenPart_grmompdgId[i] == 6:
                    top_lep = top_vec
                    top_had = anti_top_vec
                elif event.GenPart_grmompdgId[i] == -6:
                    top_lep = anti_top_vec
                    top_had = top_vec

        jets = []
        for i in range(event.nGenJetAK8):
            jets.append(TLorentzVector())
            jets[-1].SetPtEtaPhiM(event.GenJetAK8_pt[i], event.GenJetAK8_eta[i],
                                  event.GenJetAK8_phi[i], event.GenJetAK8_mass[i])
        deltas = [jets[j].DeltaR(top_had) for j in range(event.nGenJetAK8)]
        nearest_jet_idx = np.argmin(deltas)
        nearest_jet_pt = jets[nearest_jet_idx].Pt()

        jet_constituents = JetConstituents.get(event=event, index=nearest_jet_idx)

        triplets.append(Triplets.make(jet_pt=nearest_jet_pt, particle_vectors=jet_constituents))

    return np.concatenate(triplets, axis=0)


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')        # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    args = argParser.parse_args()

    all_triplets = calc_triplet_data()

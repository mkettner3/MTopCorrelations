# -*- coding: utf-8 -*-

"""
Script to determine the types of particles inside the hadronic top jet.
"""

import time
from MTopCorrelations.samples.nanoTuples_UL_RunII_nanoAOD import UL2018
from calc_triplet_data import find_hadronic_jet
from ROOT import TH1F, TLorentzVector
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse
from collections import OrderedDict
from typing import Any


def get_jet_constituents(event, level, index, max_numb_of_cons=None):
    # type: (Any, str, int, int) -> tuple
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
    cons_pdgid = []
    for i in range(event_para['Jet_cons']):
        if event_para['Jet_cons_jetIndex'][i] == index:
            particle = TLorentzVector()
            particle.SetPtEtaPhiM(event_para['Jet_cons_pt'][i], event_para['Jet_cons_eta'][i],
                                  event_para['Jet_cons_phi'][i], event_para['Jet_cons_mass'][i])
            constituents.append(particle)
            cons_pdgid.append(event_para['Jet_cons_pdgId'][i])

    def get_pt(part):
        return part.Pt()

    constituents.sort(key=get_pt, reverse=True)
    if isinstance(max_numb_of_cons, int):
        if len(constituents) > max_numb_of_cons:
            constituents = constituents[:max_numb_of_cons]

    return constituents, cons_pdgid


def calc_triplets_and_hist(samples, nbins=50, hist_range=(0, 3)):
    read_variables = [
        "nGenPart/I",
        "GenPart[pt/F,eta/F,phi/F,m/F,pdgId/I,mompdgId/I,grmompdgId/I]",
        "nGenJetAK8/I",
        "GenJetAK8[pt/F,eta/F,phi/F,mass/F]",
        "nGenJetAK8_cons/I",
        VectorTreeVariable.fromString("GenJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000),
        "nPFJetAK8/I",
        "PFJetAK8[pt/F,eta/F,phi/F,mass/F]",
        "nPFJetAK8_cons/I",
        VectorTreeVariable.fromString("PFJetAK8_cons[pt/F,eta/F,phi/F,mass/F,pdgId/I,jetIndex/I]", nMax=1000)]

    pdg_hists = [[OrderedDict() for _ in range(len(samples))] for _ in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample in enumerate(samples):
            r = sample.treeReader(variables=read_variables, selectionString="Sum$({:}JetAK8_pt>400)>=1".format(level))

            global count
            r.start()
            while r.run():                                                              # Event-Loop
                count += 1
                hadronic_jet_idx, _, _ = find_hadronic_jet(r.event, level=level, merge_tolerance=0.8, jet_pt_min=400)

                if hadronic_jet_idx is not None:
                    jet_constituents, cons_pdgId = get_jet_constituents(event=r.event, level=level,
                                                                        index=hadronic_jet_idx, max_numb_of_cons=100)
                    for pdg in cons_pdgId:
                        try:
                            pdg_hists[g][h][pdg] += 1
                        except KeyError:
                            pdg_hists[g][h][pdg] = 1

    for g, level in enumerate(['Gen', 'PF']):
        for h, sample in enumerate(samples):
            pdg_hists[g][h] = OrderedDict(sorted(pdg_hists[g][h].items(), key=lambda t: t[1], reverse=True))

    return pdg_hists


mc = [UL2018.TTbar_1, UL2018.TTbar_2, UL2018.TTbar_3]


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    args = argParser.parse_args()

    samples = [sample.split(n=args.nJobs, nSub=args.job) for sample in mc]

    start = time.time()
    count = 0
    pdg_hists = calc_triplets_and_hist(samples=samples)
    end = time.time()

    for g, level in enumerate(['Gen', 'PF']):
        print('== '+level+' ==')
        for h, sample in enumerate(samples):
            print('-- '+sample.name+' --')
            for par_type in pdg_hists[g][h].items():
                print(str(par_type[0])+': '+str(par_type[1]))
            print('')

    print('Executing determine_constituents_types.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))
    print('Number of considered events in all samples: {:}'.format(count))

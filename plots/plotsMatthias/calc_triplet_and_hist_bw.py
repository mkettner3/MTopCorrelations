# -*- coding: utf-8 -*-

"""
Script to calculate triplets and generate histograms event-wise. The triplet data is only stored in RAM for one event.
"""

import time
from math import isnan
import numpy as np
import MTopCorrelations.samples.nanoTuples_UL2018_nanoAOD as UL2018
from calc_triplet_and_hist import pt_jet_ranges
from triplet_maker import make_triplets_and_cut, make_triplets_and_cut_sim_eff
from find_hadronic_jet import find_hadronic_jet
from jet_constituents import get_jet_constituents
import ROOT
from RootTools.core.TreeVariable import VectorTreeVariable
import argparse


def cons_pt_var_func(var_fac, cons_pt):
    if cons_pt <= 10:
        return var_fac
    else:
        return 0.0375*var_fac*cons_pt + 0.625*var_fac


def calc_triplets_and_hist(sample, rew_samples, pt_jet_ranges, max_delta_zeta=float('nan'), delta_legs=float('nan'),
                           shortest_side=0.1, nbins=900, hist_range=(0, 3), pt_variations=None):
    global count, number_events_ttbar

    read_variables = [
        "Generator_weight/D",
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

    hists = [[[[ROOT.TH1D("Correlator Top ({:}, {:}, {:}, {:})".format(i, j, k, l), "3 #zeta", nbins, hist_range[0], hist_range[1]) for l in range(2)] for k in range(len(pt_jet_ranges))] for j in range(len(rew_samples))] for i in range(2)]
    hists_w = [[[[ROOT.TH1D("Correlator W ({:}, {:}, {:}, {:})".format(i, j, k, l), "3 #zeta", nbins, hist_range[0], hist_range[1]) for l in range(2)] for k in range(len(pt_jet_ranges))] for j in range(len(rew_samples))] for i in range(2)]
    hists_varied_jet = [[[[ROOT.TH1D("Correlator varied jet ({:}, {:}, {:}, {:})".format(i, j, k, l), "3 #zeta", nbins, hist_range[0], hist_range[1]) for l in range(len(pt_variations))] for k in range(len(pt_jet_ranges))] for j in range(len(rew_samples))] for i in range(2)]
    hists_varied_cons_pt = [[[[ROOT.TH1D("Correlator varied cons pt ({:}, {:}, {:}, {:})".format(i, j, k, l), "3 #zeta", nbins, hist_range[0], hist_range[1]) for l in range(len(pt_variations))] for k in range(len(pt_jet_ranges))] for j in range(len(rew_samples))] for i in range(2)]
    hists_varied_cons_eta_phi = [[[[[ROOT.TH1D("Correlator varied cons eta phi ({:}, {:}, {:}, {:}, {:})".format(i, j, k, m, n), "3 #zeta", nbins, hist_range[0], hist_range[1]) for n in range(3)] for m in range(3)] for k in range(len(pt_jet_ranges))] for j in range(len(rew_samples))] for i in range(2)]
    hists_jet_pt = [[ROOT.TH1D("Hadronic Top-Jet-p_{T} ("+str(i)+", "+str(j)+")", "Jet-p_{t}", nbins, 380, 730) for j in range(len(rew_samples))] for i in range(2)]
    hists_cons_pt = [[ROOT.TH1D("Hadronic Top-Constituents-p_{T} ("+str(i)+", "+str(j)+")", "Constituent-p_{t}", nbins, 0, 30) for j in range(len(rew_samples))] for i in range(2)]
    hists_jet_mass = [[ROOT.TH1D("Hadronic Top-Jet-mass ({:}, {:})".format(i, j), "Jet-mass", nbins, 75, 300) for j in range(len(rew_samples))] for i in range(2)]
    hists_top_mass = [[ROOT.TH1D("Hadronic Top-mass ({:}, {:})".format(i, j), "Jet-mass", nbins, 167.5, 177.5) for j in range(len(rew_samples))] for i in range(2)]
    hists_event_weight = [ROOT.TH1D("Event weights ({:})".format(i), "event-weight", nbins, -1, 1) for i in range(2)]
    hists_number_events = [ROOT.TH1I("Number of events ({:})".format(i), "Number of events", 1, 0, 1) for i in range(2)]
    hists_weighted_events = [ROOT.TH1D("Weighted number of events ({:})".format(i), "Weighted number of events", 1, 0, 1) for i in range(2)]

    for g, level in enumerate(['Gen', 'PF']):
        r = sample.treeReader(variables=read_variables, selectionString="Sum$({:}JetAK8_pt>400)>=1".format(level))
        r.start()
        while r.run():                                                              # Event-Loop
            event_weight = 60 * 831.762 * 3*0.108 * (1-3*0.108)*2 * 1000 / number_events_ttbar * r.event.Generator_weight

            hadronic_jet_idx, hadronic_jet_pt, hadronic_jet_mass, m_top_had = find_hadronic_jet(r.event, level=level,
                                                                                                merge_tolerance=0.8,
                                                                                                jet_pt_min=400)

            if hadronic_jet_idx is not None:
                jet_constituents = get_jet_constituents(event=r.event, level=level,
                                                        index=hadronic_jet_idx, particle_types='charged_par',
                                                        max_numb_of_cons=50)

                if len(jet_constituents) > 0:
                    triplets, triplets_w = construct_triplets(hadronic_jet_pt, jet_constituents,
                                                              max_delta_zeta, delta_legs, shortest_side)

                    for h, mtop_bw in enumerate(rew_samples):
                        bw_weight = calc_bw_factor(event_top_mass=m_top_had, new_top_mass=mtop_bw)

                        for k, jet_range in enumerate(pt_jet_ranges):
                            if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                for three_zeta, weight in zip(triplets[0], triplets[1]):
                                    hists[g][h][k][0].Fill(three_zeta, weight*event_weight*bw_weight)
                                    hists[g][h][k][1].Fill(three_zeta, event_weight*bw_weight)
                                for three_zeta, weight in zip(triplets_w[0], triplets_w[1]):
                                    hists_w[g][h][k][0].Fill(three_zeta, weight*event_weight*bw_weight)
                                    hists_w[g][h][k][1].Fill(three_zeta, event_weight*bw_weight)
                                break

                        hists_jet_pt[g][h].Fill(hadronic_jet_pt, event_weight*bw_weight)
                        for cons in jet_constituents:
                            hists_cons_pt[g][h].Fill(cons.Pt())
                        hists_jet_mass[g][h].Fill(hadronic_jet_mass, event_weight*bw_weight)
                        hists_top_mass[g][h].Fill(m_top_had, event_weight*bw_weight)
                    hists_event_weight[g].Fill(event_weight)
                    hists_number_events[g].Fill(0.5)
                    hists_weighted_events[g].Fill(0.5, event_weight)
                    count += 1

                    for v, var_fac in enumerate(pt_variations):
                        hadronic_jet_pt_varied = hadronic_jet_pt * var_fac

                        triplets, _ = construct_triplets(hadronic_jet_pt_varied, jet_constituents,
                                                         max_delta_zeta, delta_legs, shortest_side)

                        for h, mtop_bw in enumerate(rew_samples):
                            bw_weight = calc_bw_factor(event_top_mass=m_top_had, new_top_mass=mtop_bw)

                            for k, jet_range in enumerate(pt_jet_ranges):
                                if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                    for three_zeta, weight in zip(triplets[0], triplets[1]):
                                        hists_varied_jet[g][h][k][v].Fill(three_zeta, weight*event_weight*bw_weight)
                                    break

                    for v, var_fac in enumerate(pt_variations):
                        jet_constituents_pt_varied = [constituent * cons_pt_var_func(var_fac, constituent.Pt()) for constituent in jet_constituents]

                        triplets, _ = construct_triplets(hadronic_jet_pt, jet_constituents_pt_varied,
                                                         max_delta_zeta, delta_legs, shortest_side)

                        for h, mtop_bw in enumerate(rew_samples):
                            bw_weight = calc_bw_factor(event_top_mass=m_top_had, new_top_mass=mtop_bw)

                            for k, jet_range in enumerate(pt_jet_ranges):
                                if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                    for three_zeta, weight in zip(triplets[0], triplets[1]):
                                        hists_varied_cons_pt[g][h][k][v].Fill(three_zeta, weight*event_weight*bw_weight)
                                    break

                    for m, eff_deltaR in enumerate([0.01, 0.05, 0.1]):
                        for n, eff_probability in enumerate([2, 5, 10]):
                            triplets, _ = construct_triplets(hadronic_jet_pt, jet_constituents,
                                                             max_delta_zeta, delta_legs, shortest_side,
                                                             efficiency_simulation=True,
                                                             tracker_efficiency_deltaR=eff_deltaR,
                                                             tracker_efficiency_loss_rate=eff_probability)

                            for h, mtop_bw in enumerate(rew_samples):
                                bw_weight = calc_bw_factor(event_top_mass=m_top_had, new_top_mass=mtop_bw)

                                for k, jet_range in enumerate(pt_jet_ranges):
                                    if jet_range[0] <= hadronic_jet_pt < jet_range[1]:
                                        for three_zeta, weight in zip(triplets[0], triplets[1]):
                                            hists_varied_cons_eta_phi[g][h][k][m][n].Fill(three_zeta, weight*event_weight*bw_weight)
                                        break

    return hists, hists_w, hists_jet_pt, hists_cons_pt, hists_jet_mass, hists_top_mass, hists_varied_jet, hists_varied_cons_pt, hists_varied_cons_eta_phi, hists_event_weight, hists_number_events, hists_weighted_events


def construct_triplets(hadronic_jet_pt, jet_constituents, max_delta_zeta, delta_legs, shortest_side,
                       efficiency_simulation=False, tracker_efficiency_deltaR=None, tracker_efficiency_loss_rate=None):
    # type: (float, list, float, float, float, bool, float, int) -> tuple

    if isnan(max_delta_zeta):
        max_delta_zeta_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
    else:
        max_delta_zeta_calc = max_delta_zeta

    if isnan(delta_legs):
        delta_legs_calc = 3.5 / 3. * (170. / hadronic_jet_pt) ** 2
    else:
        delta_legs_calc = delta_legs

    if not efficiency_simulation:
        (triplets,
         triplets_w) = make_triplets_and_cut(jet_pt=hadronic_jet_pt, particle_vectors=jet_constituents,
                                             max_delta_zeta=max_delta_zeta_calc,
                                             delta_legs=delta_legs_calc, shortest_side=shortest_side)
    else:
        (triplets,
         triplets_w) = make_triplets_and_cut_sim_eff(jet_pt=hadronic_jet_pt, particle_vectors=jet_constituents,
                                                     max_delta_zeta=max_delta_zeta_calc,
                                                     delta_legs=delta_legs_calc, shortest_side=shortest_side,
                                                     tracker_efficiency_deltaR=tracker_efficiency_deltaR,
                                                     tracker_efficiency_loss_rate=tracker_efficiency_loss_rate)

    return triplets, triplets_w


def calc_bw_factor(event_top_mass, new_top_mass):
    if new_top_mass is not None:
        factor = breit_wigner(1.324, new_top_mass, event_top_mass) / breit_wigner(1.324, 172.5, event_top_mass)
    else:
        factor = 1

    return factor


def breit_wigner(width, peak_position, mass):
    gamma = np.sqrt(peak_position**2 * (peak_position**2 + width**2))
    k = (2*np.sqrt(2)*peak_position*width*gamma) / (np.pi*np.sqrt(peak_position**2 + gamma))
    value = k / ((mass**2-peak_position**2)**2 + peak_position**2 * width**2)

    return value


def save_root_hists(hists_top, hists_w, hists_jet_pt, hists_cons_pt, hists_jet_mass, hists_top_mass, hists_varied_jet,
                    hists_varied_cons_pt, hists_varied_cons_eta_phi, pt_variations, hists_ev_weight,
                    hists_number_events, hists_weighted_events, rew_values, pt_jet_ranges, filename):
    # type: (list, list, list, list, list, list, list, list, list, list, list, list, list, list, list, str) -> None

    r_file = ROOT.TFile(filename, 'RECREATE')
    r_file.cd()
    for p_type in ['Top-Quark', 'W-Boson', 'Others']:
        for lev in ['Gen-Level', 'PF-Level']:
            if p_type in ['Top-Quark', 'W-Boson']:
                for weight in ['weighted', 'absolute']:
                    r_file.mkdir(p_type+'/'+lev+'/'+weight)
            else:
                r_file.mkdir(p_type+'/'+lev)

    for g, level in enumerate(['Gen', 'PF']):
        for h, mtop_bw in enumerate(rew_values):
            for k, pt_jet_range in enumerate(pt_jet_ranges):
                for f, (weighted, weight_dir) in enumerate(zip(['', '_abscou'], ['weighted', 'absolute'])):
                    r_file.cd('/Top-Quark/'+level+'-Level/'+weight_dir)
                    hists_top[g][h][k][f].Write('correlator_hist_{:}_{:}_{:}_{:}{:}'.format(level, mtop_bw, pt_jet_range[0], pt_jet_range[1], weighted))
                    r_file.cd('/W-Boson/'+level+'-Level/'+weight_dir)
                    hists_w[g][h][k][f].Write('correlator_hist_W_{:}_{:}_{:}_{:}{:}'.format(level, mtop_bw, pt_jet_range[0], pt_jet_range[1], weighted))
                for v, var_fac in enumerate(pt_variations):
                    r_file.cd('/Top-Quark/'+level+'-Level/weighted')
                    hists_varied_jet[g][h][k][v].Write('correlator_hist_varied_jet_{:.2f}_{:}_{:}_{:}_{:}'.format(var_fac, level, mtop_bw, pt_jet_range[0], pt_jet_range[1]))
                    hists_varied_cons_pt[g][h][k][v].Write('correlator_hist_varied_cons_pt_{:.2f}_{:}_{:}_{:}_{:}'.format(var_fac, level, mtop_bw, pt_jet_range[0], pt_jet_range[1]))
                for m, eff_deltaR in enumerate([0.01, 0.05, 0.1]):
                    for n, eff_probability in enumerate([2, 5, 10]):
                        hists_varied_cons_eta_phi[g][h][k][m][n].Write('correlator_hist_varied_cons_eta_phi_{:}_{:}_{:}_{:}_{:}_{:}'.format(level, mtop_bw, pt_jet_range[0], pt_jet_range[1], eff_deltaR, eff_probability))
        r_file.cd('/Others/'+level+'-Level')
        for h, mtop_bw in enumerate(rew_values):
            hists_jet_pt[g][h].Write('hadronic_top_jet_pt_hist_{:}_{:}'.format(level, mtop_bw))
            hists_cons_pt[g][h].Write('hadronic_top_constituents_pt_hist_{:}_{:}'.format(level, mtop_bw))
            hists_jet_mass[g][h].Write('hadronic_top_jet_mass_hist_{:}_{:}'.format(level, mtop_bw))
            hists_top_mass[g][h].Write('hadronic_top_mass_hist_{:}_{:}'.format(level, mtop_bw))
        hists_ev_weight[g].Write('event_weights_{:}'.format(level))
        hists_number_events[g].Write('number_of_events_{:}'.format(level))
        hists_weighted_events[g].Write('weighted_number_of_events_{:}'.format(level))

    r_file.Close()


sample = UL2018.TTbar_3
number_events_ttbar = 139438020309

rew_samples = [171.5, 171.75, 172., 172.25, None, 172.75, 173., 173.25, 173.5]   # Breit-Wigner reweighted samples


if __name__ == '__main__':
    argParser = argparse.ArgumentParser(description='Argument parser')    # So #SPLIT100 can be used in the bash script.
    argParser.add_argument('--nJobs', action='store', nargs='?', type=int, default=1)
    argParser.add_argument('--job', action='store', type=int, default=0)
    argParser.add_argument('--pt_variations', action='store', type=str, default='1.1, 1.05, 1.02, 1.01, 0.99, 0.98, 0.95, 0.9')
    args = argParser.parse_args()

    sample = sample.split(n=args.nJobs, nSub=args.job)
    pt_variations = [float(e.strip()) for e in args.pt_variations.split(',')]

    start = time.time()
    count = 0
    (hists, hists_w,
     hists_jet_pt, hists_cons_pt, hists_jet_mass, hists_top_mass,
     hists_varied_jet, hists_varied_cons_pt, hists_varied_cons_eta_phi,
     hists_event_weight, hists_number_events, hists_weighted_events) = calc_triplets_and_hist(
                                                  sample=sample, rew_samples=rew_samples, pt_jet_ranges=pt_jet_ranges,
                                                  max_delta_zeta=0.13, delta_legs=float('nan'), shortest_side=0.05,
                                                  pt_variations=pt_variations)
    save_root_hists(hists_top=hists, hists_w=hists_w, hists_jet_pt=hists_jet_pt, hists_cons_pt=hists_cons_pt,
                    hists_jet_mass=hists_jet_mass, hists_top_mass=hists_top_mass, hists_varied_jet=hists_varied_jet,
                    hists_varied_cons_pt=hists_varied_cons_pt, hists_varied_cons_eta_phi=hists_varied_cons_eta_phi,
                    pt_variations=pt_variations,
                    hists_ev_weight=hists_event_weight, hists_number_events=hists_number_events,
                    hists_weighted_events=hists_weighted_events, rew_values=rew_samples, pt_jet_ranges=pt_jet_ranges,
                    filename='histogram_files/correlator_hist_trip_27_pp_{:03}.root'.format(args.job))
    end = time.time()

    print('Executing calc_triplet_and_hist.py took {:.0f}:{:.2f} min:sec.'.format((end-start)//60, (end-start)%60))
    print('Number of considered events in all samples: {:}'.format(count))

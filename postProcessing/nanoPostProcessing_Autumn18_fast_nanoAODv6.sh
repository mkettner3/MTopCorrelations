
# 1lep (gamma) -> flag ttbar/ttgamma train fisher information in 1l final state with and w/o background

# 2l ttZ / DY+bb 

#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample ttG_noFullyHad_fast #SPLIT500
#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample WGToLNu_fast #SPLIT200
#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample ZGTo2L_fast #SPLIT200
#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample ttW01j_fast #SPLIT100
#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample WZTojj2L_fast #SPLIT 99
#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample WZToLNujj_fast #SPLIT100
#python nanoPostProcessing.py  --overwrite --forceProxy --skim singlelep-photon --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WW_fast #SPLIT100


python nanoPostProcessing.py  --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample ttZ01j_fast #SPLIT100
python nanoPostProcessing.py  --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample WZTo3L1Nu_fast #SPLIT100
python nanoPostProcessing.py  --overwrite --forceProxy --skim trilep --year 2018 --processingEra tWZ_nAODv6_private_v6 --sample ZZ_fast #SPLIT100

#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample ttZ01j_fast #SPLIT98
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZTo3L1Nu_fast #SPLIT100
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZTojj2L_fast #SPLIT 99
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZToLNujj_fast #SPLIT100
#
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_LO #SPLIT40
#
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --LHEHTCut 70 --sample DYJetsToLL_M50_LO #SPLIT100
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT70to100 #SPLIT50
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT100to200 #SPLIT50
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT200to400 #SPLI50
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT400to600 #SPLIT100
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT600to800 #SPLIT50
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT800to1200 #SPLIT30
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT10
#python nanoPostProcessing.py  --overwrite --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT2500toInf #SPLIT10

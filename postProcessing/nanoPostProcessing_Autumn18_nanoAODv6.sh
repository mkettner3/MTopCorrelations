## DY

## full stats
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_LO #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M10to50_LO #SPLIT29
#
## HT binned samples ##
## high mass #
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --LHEHTCut 70 --sample DYJetsToLL_M50_LO #SPLIT40
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT70to100 #SPLIT16
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT100to200 #SPLIT12
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT200to400 #SPLIT12
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT400to600 #SPLIT20
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT600to800 #SPLIT13
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT800to1200 #SPLIT6
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT1200to2500 #SPLIT2
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M50_HT2500toInf #SPLIT2
## low mass #
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --LHEHTCut 70 --sample DYJetsToLL_M10to50_LO #SPLIT28
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M4to50_HT70to100 #SPLIT16
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M4to50_HT100to200 #SPLIT11
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M4to50_HT200to400 #SPLIT2
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M4to50_HT400to600 #SPLIT3
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample DYJetsToLL_M4to50_HT600toInf #SPLIT5

# top
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --flagTTBar --sample TTLep_pow  #SPLIT250
##python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --flagTTBar --sample TTSingleLep_pow #SPLIT80
##python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --flagTTGamma --sample TTGLep #SPLIT20
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --reduceSizeBy 5 --sample TToLeptons_sch_amcatnlo #SPLIT5
#python nanoPostProcessing.py --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample T_tWch #SPLIT3
#python nanoPostProcessing.py --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TBar_tWch #SPLIT3
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --reduceSizeBy 15 --sample T_tch_pow #SPLIT8
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --reduceSizeBy 15 --sample TBar_tch_pow #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tZq_ll #SPLIT17
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll 
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTWToLNu #SPLIT50
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTWToQQ #SPLIT10
###python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTW_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTZToLLNuNu #SPLIT100
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTZToLLNuNu_m1to10 #SPLIT5
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTZToQQ #SPLIT16
###python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTZ_LO #SPLIT20
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTHbbLep #SPLIT100
python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTHnobb_pow #SPLIT120
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample THQ #SPLIT4
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample THW #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWnunu 

## di/multi boson

#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample VVTo2L2Nu #SPLIT13
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZTo3LNu_amcatnlo #SPLIT16
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZTo2L2Q #SPLIT20
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZTo1L3Nu #SPLIT2
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WWTo1L1Nu2Q #SPLIT8
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WWToLNuQQ #SPLIT18
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample ZZTo4L #SPLIT8
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample ZZTo2Q2Nu #SPLIT10
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --reduceSizeBy 5 --sample ZZTo2L2Q #SPLIT5
## additional samples for ARC studies
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WWTo2L2Nu #SPLIT11
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample ZZTo2L2Nu #SPLIT11
###
####python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WW #SPLIT5
####python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample ZZ #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WWW_4F #SPLIT3
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WWZ 
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample WZZ #SPLIT2
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample ZZZ 
##
### rare
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTTT #SPLIT5
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTWZ 
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTWW #SPLIT3
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample TTZZ 
#
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll_thad_Wlept_DR #SPLIT10
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll_thad_Wlept_DS #SPLIT10
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll_tlept_Whad_DR #SPLIT10
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll_tlept_Whad_DS #SPLIT10
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll_tlept_Wlept_DR #SPLIT10
#python nanoPostProcessing.py  --forceProxy --skim dilep --year 2018 --processingEra tWZ_nAODv6_private_v4 --sample tWll_tlept_Wlept_DS #SPLIT10

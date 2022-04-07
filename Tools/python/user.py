import os


if os.environ["USER"] in ["dennis.schwarz"]:
    postprocessing_output_directory = "/scratch-cbe/users/dennis.schwarz/MTopCorrelations/nanoTuples"
    postprocessing_tmp_directory    = "/scratch/hephy/cms/dennis.schwarz/MTopCorrelations/tmp/"
    plot_directory                  = "/groups/hephy/cms/dennis.schwarz/www/MTopCorrelations/plots"
    cache_dir                       = "/groups/hephy/cms/dennis.schwarz/MTopCorrelations/caches"
    analysis_results                = "/groups/hephy/cms/dennis.schwarz/MTopCorrelations/results/v1"
    mva_directory                   = "/groups/hephy/cms/dennis.schwarz/MTopCorrelations/MVA"
    cern_proxy_certificate          = "/users/dennis.schwarz/.private/.proxy"
    combineReleaseLocation          = "/users/dennis.schwarz/CMSSW_10_6_0/src"

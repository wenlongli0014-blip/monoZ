
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do compute_instrMET_weights.py -o data/InstrMetReweighting/weight_nvtx_${YEAR}.root -s nvtx batch_dilepton_${YEAR}/merged/ batch_singlephoton_${YEAR}/merged/; done

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_nvtx_reweight --config singlephoton/${YEAR}-ul-nvtxReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_nvtx_reweight && bash send_jobs.sh); done

sleep 10m

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_singlephoton_${YEAR}_nvtx_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt; done

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do compute_instrMET_weights.py -o data/InstrMetReweighting/weight_eta_${YEAR}.root -s eta batch_dilepton_${YEAR}/merged/ batch_singlephoton_${YEAR}_nvtx_reweight/merged/; done

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_eta_reweight --config singlephoton/${YEAR}-ul-etaReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_eta_reweight && bash send_jobs.sh); done

sleep 10m

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_singlephoton_${YEAR}_eta_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt; done

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do compute_instrMET_weights.py -o data/InstrMetReweighting/weight_pt_${YEAR}.root -s pt batch_dilepton_${YEAR}/merged/ batch_singlephoton_${YEAR}_eta_reweight/merged/; done

# compute_mass_lineshape.py -o data/InstrMetReweighting/lineshape_mass_2018.root  batch_dilepton_2018/merged/

#计算Zgamma的normalization并reweight
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python ZGammaReight.py -y ${YEAR}; done

#进入另一个环境
. ./env_reweight.sh

#计算WJets的fakeRatio

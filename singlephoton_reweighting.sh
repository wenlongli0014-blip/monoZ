
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_nvtx_2018.root -s nvtx batch_dilepton_2018/merged/ batch_singlephoton_2018/merged/

for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_nvtx_reweight --config singlephoton/${YEAR}-ul-nvtxReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_nvtx_reweight && bash send_jobs.sh); done
sleep 360
harvest.py --task-dir batch_singlephoton_2018_nvtx_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018_data.txt

compute_instrMET_weights.py -o data/InstrMetReweighting/weight_eta_2018.root -s eta batch_dilepton_2018/merged/ batch_singlephoton_2018_nvtx_reweight/merged/
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_eta_reweight --config singlephoton/${YEAR}-ul-etaReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_eta_reweight && bash send_jobs.sh); done
sleep 360
harvest.py --task-dir batch_singlephoton_2018_eta_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018_data.txt

compute_instrMET_weights.py -o data/InstrMetReweighting/weight_pt_2018.root -s pt batch_dilepton_2018/merged/ batch_singlephoton_2018_eta_reweight/merged/

compute_mass_lineshape.py -o data/InstrMetReweighting/lineshape_mass_2018.root  batch_dilepton_2018/merged/
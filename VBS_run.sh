#test
runHZZanalysis --config singlephoton/2018-ul.yaml --ddf /pnfs/iihe/cms/store/user/yunyangl/ULsamples/photon/2023-08-22_2018-NanoAODv9/DDF/InstrMET/QCD_HT50to100.yaml --analysis PhotonTrees --more-vars
#dilepton,singlephoton,egammafrommisid,singleelectron,zgamma
# process samples
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_dilepton_${YEAR} --config dilepton/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_dilepton_${YEAR}.txt --analysis DileptonTrees --more-vars && (cd batch_dilepton_${YEAR} && bash send_jobs.sh); done
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR} --config singlephoton/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoto_${YEAR} && bash send_jobs.sh); done
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_egammafrommisid_${YEAR} --config egammafrommisid/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_EGammaFromMisid_${YEAR}.txt --analysis EGammaFromMisid --more-vars && (cd batch_egammafrommisid_${YEAR} && bash send_jobs.sh); done
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singleelectron_${YEAR} --config singleelectron/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_SingleElectron_${YEAR}.txt --analysis ElectronTrees --more-vars && (cd batch_singleelectron_${YEAR} && bash send_jobs.sh); done
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_zgamma_${YEAR} --config zgamma/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_ZGamma_dilepton_${YEAR}.txt --analysis ZGammaTrees --more-vars && (cd batch_zgamma_${YEAR} && bash send_jobs.sh); done

#merge the rootfile 
harvest.py --task-dir batch_dilepton_2018/ --config 2018-ul.yaml  $HZZ2L2NU_BASE/config/samples_dilepton_2018.txt
harvest.py --task-dir batch_singlephoton_2018/ --config 2018-ul.yaml  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018.txt
harvest.py --task-dir batch_egammafrommisid_2018/ --config 2018-ul.yaml  $HZZ2L2NU_BASE/config/samples_EGammaFromMisid_2018.txt
harvest.py --task-dir batch_singleelectron_2018/ --config 2018-ul.yaml  $HZZ2L2NU_BASE/config/samples_SingleElectron_2018.txt
harvest.py --task-dir batch_zgamma_2018/ --config 2018-ul.yaml  $HZZ2L2NU_BASE/config/samples_ZGamma_dilepton_2018.txt

#plot
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_dilepton.yaml --prefix ./batch_dilepton_2018/merged/ --output data_sim_dilepton_2018 --year 2018
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_photon.yaml --prefix ./batch_singlephoton_2018/merged/ --output data_sim_singlephoton_2018 --year 2018
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_egammafrommisid.yaml --prefix ./batch_egammafrommisid_2018/merged/ --output data_sim_egammafrommisid_2018 --year 2018
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_electron.yaml --prefix ./batch_singleelectron_2018/merged/ --output data_sim_singleelectron_2018 --year 2018
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_zgamma.yaml --prefix ./batch_zgamma_2018/merged/ --output data_sim_zgamma_2018 --year 2018

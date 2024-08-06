#test
# runHZZanalysis --config singlephoton/2018-ul.yaml --ddf /pnfs/iihe/cms/store/user/yunyangl/ULsamples/photon/2023-08-22_2018-NanoAODv9/DDF/InstrMET/QCD_HT50to100.yaml --analysis PhotonTrees --more-vars

# process samples

for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_dilepton_${YEAR} --config dilepton/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_dilepton_${YEAR}.txt --analysis DileptonTrees --more-vars && (cd batch_dilepton_${YEAR} && bash send_jobs.sh); done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR} --config singlephoton/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR} && bash send_jobs.sh); done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_egammafrommisid_${YEAR} --config egammafrommisid/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_EGammaFromMisid_${YEAR}.txt --analysis EGammaFromMisid --more-vars && (cd batch_egammafrommisid_${YEAR} && bash send_jobs.sh); done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_singleelectron_${YEAR} --config singleelectron/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_SingleElectron_${YEAR}.txt --analysis ElectronTrees --more-vars && (cd batch_singleelectron_${YEAR} && bash send_jobs.sh); done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do prepare_htcondor_jobs.py --task-dir batch_zgamma_${YEAR} --config zgamma/${YEAR}-ul.yaml -- $HZZ2L2NU_BASE/config/samples_ZGamma_dilepton_${YEAR}.txt --analysis ZGammaTrees --more-vars && (cd batch_zgamma_${YEAR} && bash send_jobs.sh); done

sleep 4h
#merge the rootfile 
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_dilepton_${YEAR}/ --config ${YEAR}-ul.yaml  $HZZ2L2NU_BASE/config/samples_dilepton_${YEAR}.txt; done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_singlephoton_${YEAR}/ --config ${YEAR}-ul.yaml  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}.txt; done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_egammafrommisid_${YEAR}/ --config ${YEAR}-ul.yaml  $HZZ2L2NU_BASE/config/samples_EGammaFromMisid_${YEAR}.txt; done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_singleelectron_${YEAR}/ --config ${YEAR}-ul.yaml  $HZZ2L2NU_BASE/config/samples_SingleElectron_${YEAR}.txt; done
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do harvest.py --task-dir batch_zgamma_${YEAR}/ --config ${YEAR}-ul.yaml  $HZZ2L2NU_BASE/config/samples_ZGamma_dilepton_${YEAR}.txt; done

#plot
# for YEAR in 2016HIPM 2016noHIPM 2017 2018; do plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_dilepton.yaml --prefix ./batch_dilepton_${YEAR}/merged/ --output data_sim_dilepton_${YEAR} --year 2018; done
# for YEAR in 2016HIPM 2016noHIPM 2017 2018; do plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_photon.yaml --prefix ./batch_singlephoton_${YEAR}/merged/ --output data_sim_singlephoton_${YEAR} --year 2018; done
# for YEAR in 2016HIPM 2016noHIPM 2017 2018; do plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_egammafrommisid.yaml --prefix ./batch_egammafrommisid_${YEAR}/merged/ --output data_sim_egammafrommisid_${YEAR} --year 2018; done
# for YEAR in 2016HIPM 2016noHIPM 2017 2018; do plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_electron.yaml --prefix ./batch_singleelectron_${YEAR}/merged/ --output data_sim_singleelectron_${YEAR} --year 2018; done
# for YEAR in 2016HIPM 2016noHIPM 2017 2018; do plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_zgamma.yaml --prefix ./batch_zgamma_${YEAR}/merged/ --output data_sim_zgamma_${YEAR} --year 2018; done

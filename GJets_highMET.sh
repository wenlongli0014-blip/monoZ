#退出当前会话，重新打开，进入另一个环境
. ./env_reweight.sh

hadd batch_singlephoton_2016/merged/ZNuNuGJets.root batch_singlephoton_2016HIPM/merged/ZNuNuGJets.root batch_singlephoton_2016noHIPM/merged/ZNuNuGJets.root

cd zz2l2nu_reweight

#计算Zgamma的normalization并reweight
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python ZGammaReweight.py -y ${YEAR}; done

#计算WJets的fakeRatio
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python fakeratio.py -y ${YEAR}; done

#计算reweight后的WJets
cd WJets
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python WJets.py -y ${YEAR}; done
cd ..

#计算扣除本底后的GJets
cd InstrMETReweighting
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python CollectWeights.py -y ${YEAR}; done
cd ..
#计算Data Driven DY在low MET 和 SR 的估计值
python Data_Driven_DY_2018.py
for YEAR in 2016HIPM 2016noHIPM 2017; do python Data_Driven_DY.py -y ${YEAR}; done
cd root
hadd histo2d_2016_SR.root histo2d_2016HIPM_SR.root histo2d_2016noHIPM_SR.root 
hadd histo2d_2016_lowMET.root histo2d_2016HIPM_lowMET.root histo2d_2016noHIPM_lowMET.root   
hadd histo2d_2017_and_2018_lowMET.root histo2d_2017_lowMET.root histo2d_2018_lowMET.root 
cd ../..
mkdir batch_singlephoton_2016
mkdir batch_singlephoton_2016/merged
hadd batch_singlephoton_2016/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root  batch_singlephoton_2016HIPM/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root batch_singlephoton_2016noHIPM/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root
hadd batch_singlephoton_2017_and_2018/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root  batch_singlephoton_2017/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root batch_singlephoton_2018/merged/DYJetsToLL_Data-driven_GJets_nvtx_eta_pt_reweighted.root
cd zz2l2nu_reweight

#将pt_ll[60,82.5] MET>160,200区间的DY事例数替换成MC
for YEAR in 2016 2017; do python replaceGJets_2016_2017.py -y ${YEAR}; done
python replaceGJets_2018.py 
hadd histo2d_2017_and_2018_SR_final.root histo2d_2017_SR_final.root histo2d_2018_SR_final.root
python final_plot.py
#退出当前会话，重新打开，进入另一个环境
. ./env_reweight.sh
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
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python InstrMETReweighting/CollectWeights.py -y ${YEAR}; done
cd ..
#计算Data Driven DY在low MET 和 SR 的估计值
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python Data_Driven_DY.py -y ${YEAR}; done

#将pt_ll[60,82.5] MET>160区间的DY事例数替换成MC
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python replaceGJets.py -y ${YEAR}; done
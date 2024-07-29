#退出当前会话，重新打开，进入另一个环境
. ./env_reweight.sh
cd zz2l2nu_reight
#计算Zgamma的normalization并reweight
for YEAR in 2016HIPM 2016noHIPM 2017 2018; do python ZGammaReight.py -y ${YEAR}; done

#计算WJets的fakeRatio
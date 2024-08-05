# Consistent LCG environment (http://lcginfo.cern.ch)
. /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_105 x86_64-centos7-gcc11-opt
export SCRAM_ARCH=slc7_amd64_gcc920

export HZZ2L2NU_BASE=$(pwd)
export PYTHONPATH="${HZZ2L2NU_BASE}/python:$PYTHONPATH"
export PATH="${HZZ2L2NU_BASE}/bin:${PATH}"

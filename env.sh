source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_105 x86_64-el9-gcc12-opt
export SCRAM_ARCH=el9_amd64_gcc12

export HZZ2L2NU_BASE=$(pwd)
export PYTHONPATH="${HZZ2L2NU_BASE}/python:$PYTHONPATH"
export PATH="${HZZ2L2NU_BASE}/bin:${PATH}"

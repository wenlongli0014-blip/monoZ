## What `HZZ` framework do?
This framework can be used to process NanoAOD samples, including processing various variables, correction, and selection.

## Set up and building
This repository cloned with
```
git clone ssh://git@gitlab.cern.ch:7999/hgao/hzz2l2nu.git .
```

At the start of each session, set up the environment with
```sh
. ./env.sh
```

This script also stores the path to the base directory in environment variable `HZZ2L2NU_BASE`, which should then be used in scripts and compiled code to resolve relative paths to auxiliary data files, such as data-driven weights.

Build the package with the following commands:

```sh
rm -rf python/hzz/pyroothist/
git clone https://github.com/andrey-popov/pyroothist.git python/hzz/pyroothist
mkdir build
cd build
cmake ..
make -j $(nproc)
```

The warning from CMake about the new version of Boost can be safely ignored. Executable `runHZZanalysis` is put into `$HZZ2L2NU_BASE/bin`, and it is accessible from `$PATH`. To rebuild the package after a change has been introduced to the code, repeat `make`. To start the build from scratch, remove the directory `build` and repeat the commands above.

It is also possible to create a program outside of the repository and link it against the shared library of the framework. See [here](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/wikis/shared-library) for documentation.


1. 5 regions including `Dilepton, EGammaFromMisid, SingleElectron, SinglePhoton, and ZGamma` have been designed for signal selection and estimation of DY Data Driven in the VBS ZZ analysis.
- config/samples_*.txt provides samples.
- config/REGION(eg.dilepton)/*.yaml provides triggers,corrections and selections.
- src/*Trees.cc defines the variables to be stored.
- You can design new regions based on your own analysis by adding the above files.
2. 
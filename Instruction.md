## What `HZZ` framework do?
This framework can be used to process NanoAOD samples, including processing various variables, correction, selection and plotting variable distribution.

5 regions including `Dilepton, EGammaFromMisid, SingleElectron, SinglePhoton, and ZGamma` have been designed for signal selection and estimation of DY Data Driven in the VBS ZZ analysis.
- config/samples_*.txt provides samples.
- config/REGION(eg.dilepton)/*.yaml provides triggers,corrections and selections.
- src/*Trees.cc defines the variables to be stored.
- You can design new regions based on your own analysis by adding the above files.
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

The warning from CMake about the new version of Boost can be safely ignored. Executable `runHZZanalysis` is put into `$HZZ2L2NU_BASE/bin`, and it is accessible from `$PATH`.

To rebuild the package after a change has been introduced to the code, repeat the commands
```
(cd build && make -j $(nproc))
```
## Running in local
Computationally heavy part of the analysis is carried out by program runHZZanalysis. Here is an example command to run it interactively:
```
runHZZanalysis --config dilepton/2018-ul.yaml \
--ddf /pnfs/iihe/cms/store/user/hanwen/DileptonUL/2023-09-11_2016HIPM-NanoAODv9/DDF/Dilepton/DYJetsToLL_PtZ-0To50.yaml
--analysis DileptonTrees \
--max-events 10000 \
--more-vars \
```
- The first parameter `config` is the path to the master configuration file, such as 2016.yaml. It provides global settings that affect all analyises and all datasets. The path is resolved with the help of FileInPath service. Standard configuration files are located in directory $HZZ2L2NU_BASE/config, which is checked by FileInPath automatically.
- The second parameter `ddf` is the path to a dataset definition file (either a full one or a derived fragment). It provides paths to input files included in the dataset and all dataset-specific configuration parameters.
- The third parameter `analysis` specify which analysis should be executed.
- The fourth parameter `max-events` is the maximal number of events to process.
- The fifth parameter `more-vars`is to add more variables to the processed root.

A number of other command line parameters are supported, many of them also have shortcuts. The complete list can be obtained by running
```
runHZZanalysis --analysis <analysis> --help
```
## Running in local

# What `HZZ` framework do?
This framework can be used to process NanoAOD samples, including processing various variables, correction, selection and plotting variable distribution.

5 regions including `Dilepton, EGammaFromMisid, SingleElectron, SinglePhoton, and ZGamma` have been designed for signal selection and estimation of DY Data Driven in the VBS ZZ analysis.
- config/samples_*.txt provides samples.
- config/REGION(eg.dilepton)/*.yaml provides triggers,corrections and selections.
- src/*Trees.cc defines the variables to be stored.
- You can design new regions based on your own analysis by adding the above files.
## Set up and building
This repository cloned with
```
git clone ssh://git@gitlab.cern.ch:7999/hgao/hzz2l2nu.git 
```

At the start of each session, set up the environment with (you can comment out different [lines](env.sh) based on the work cluster)
```sh
cd hzz2l2nu
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
--ddf /eos/cms/store/group/phys_smp/ZZTo2L2Nu/HZZsample/2018/YAML/DYJetsToLL_PtZ-0To50.yaml \
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
## Using batch system
This repository includes scripts to submit jobs to the condor and to merge their output. [Eamples](VBS_run.sh) that include all process. 
### Submission
The submission is done with commands like 
```
prepare_htcondor_jobs.py --task-dir batch_dilepton_2018 \
--config dilepton/2018-ul.yaml \
-- $HZZ2L2NU_BASE/config/samples_dilepton_2018.txt \
--analysis DileptonTrees && (cd batch_dilepton_2018 && bash send_jobs.sh)
```
- The first parameter is the name of the directory in which job scripts and output files will be stored. 
- As with the interactive running, the second parameter is the master configuration file.
- It is also automatically forwarded to runHZZanalysis. The first positional argument is a list of datasets to be processed. Each job will process at maximum one dataset, and large datasets are automatically split among multiple jobs. All subsequent positional arguments are forwarded to the executable without a change. When some of them start with dashes, as in the example above, it is convenient to denote the start of positional arguments with `--` so that they are not interpreted as optional arguments for `prepare_htcondor_jobs.py`.
By default, `prepare_htcondor_jobs.py` will run on nodes of the batch system executable runHZZanalysis. But it is possible to supply a different executable with option `--prog`.

## Systematic variations
The full list of options supported by prepare_jobs.py is available through its help message. One of them is option `--syst` that requests systematic variations to be applied. All supported systematic uncertainties are listed in file [config/syst.yaml](https://gitlab.cern.ch/hgao/hzz2l2nu/-/blob/vbs_cuts/config/syst.yaml). For each uncertainty it provides a sequence of masks that defines which datasets are affected by this uncertainty. The masks are checked against the names of the datasets, as specified in [dataset definition files](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/wikis/dataset-definitions). Currently, there is only one way that is work:
- `--syst jec_up` or `--syst jec_down`  These are examples of fully specified variations. They will be propagated to the underlying executable without a change, but only datasets that are affected by those variations (as indicated in [config/syst.yaml](https://gitlab.cern.ch/hgao/hzz2l2nu/-/blob/vbs_cuts/config/syst.yaml)) will be processed.
- You need specific sources of systematic and shift direction.

## Harvesting
When all jobs have finished (which can be checked with `condor_q`), their outputs can be merged with
```
harvest.py --task-dir batch_dilepton_2018/ --config 2018-ul.yaml  $HZZ2L2NU_BASE/config/samples_dilepton_2018.txt
```
The will be placed in `batch_dilepton_2018/merged`. You will probably want to move that directory somewhere and delete the rest of directory `batch_dilepton_2018`.

## Post-processing with event-based analysis

Instructions below are for the event-based analyses, such as `DileptonTrees`, only. They assume that the harvesting has finished, and the merged files are available in directory `$tree_dir`.

### Plots with comparison between data and simulation

These plots can be produced by running

```sh
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim_dilepton.yaml \
--prefix ./batch_dilepton_2018/merged/ \
--output data_sim_dilepton_2018 --year 2018
```


- Event selection and variables to be plotted are described in the configuration file given as the first argument to the script. 
- The [file](config/plot_data_sim_dilepton.yaml) included in the repository is considered an example, which you would adjust to your needs. 
- Names (or rather name fragments) of input ROOT files are specified in the configuraiton. Each name is _textually_ prepended with the prefix given with the corresponding command line option  (or, alternatively, in the configuration file), which allows to specify the directory containing the files but also a common starting name fragment. 
- In addition to this prefix, multiple standard endings are tried for each file name. The directory in which produced figures will be stored is given by flag `--output`.

# DY Data Driven Estimation 
## Weights for the photon control region
[Eamples](singlephoton_reweighting.sh) that compute weight in singlephoton.
The following procedure needs to be used to re-compute the weights for the photon CR (in order to get an estimate of the Z+jets contribution in the signal region):

* Run once the `DileptonTrees` analysis (this can be done on data only), using the option `--ptmiss-cut=0`. Harvest the output.
* Run once the `PhotonTrees` analysis (this can be done on data only), changing the config file for NOT applying mean weights (else, it would drop events with a mean weight of 0).
以上两步已经在之前完成。
1.  Compute nvtx weights:

```sh
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_nvtx_2018.root -s nvtx batch_dilepton_2018/merged/ batch_singlephoton_2018/merged/
```


* Re-run the `PhotonTrees` analysis (the same way as before), making sure that nvtx weights are applied (check that in the config file).
```sh
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_nvtx_reweight --config singlephoton/${YEAR}-ul-nvtxReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_nvtx_reweight && bash send_jobs.sh); done
harvest.py --task-dir batch_singlephoton_2018_nvtx_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018_data.txt
```
2. Compute eta weights, the same way as before:
```sh
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_eta_2018.root -s eta batch_dilepton_2018/merged/ batch_singlephoton_2018_nvtx_reweight/merged/
```

* Re-run the `PhotonTrees` analysis once again, applying the nvtx and eta weights.
```sh
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_eta_reweight --config singlephoton/${YEAR}-ul-etaReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_eta_reweight && bash send_jobs.sh); done
harvest.py --task-dir batch_singlephoton_2018_eta_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018_data.txt
```
3. Compute pT weights. Notice that pT weights also include normalization:

```sh
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_pt_2018.root -s pt batch_dilepton_2018/merged/ batch_singlephoton_2018_eta_reweight/merged/
```






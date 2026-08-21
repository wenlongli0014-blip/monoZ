# mono-Z / H&rarr;ZZ&rarr;2&ell;2&nu; analysis

This repository preserves the basic structure of the upstream
[`yujil/hzz2l2nu`](https://gitlab.cern.ch/yujil/hzz2l2nu) framework and adds
the mono-Z analysis selections, signal configurations, DNN workflow, and CMS
Combine statistical workflow used in this analysis.

The GitHub and CERN GitLab repositories are mirrors of the same source state.
This repository uses Git submodules and should be cloned with one of

```sh
git clone --recurse-submodules https://github.com/wenlongli0014-blip/monoZ.git
# CERN mirror and development branch:
git clone --branch offshellH --recurse-submodules \
  https://gitlab.cern.ch/liwe/hzz2l2nu.git monoZ
```

## mono-Z extensions

Relative to the upstream framework, this repository includes:

- corrected MET covariance handling and corrected missing-momentum
  significance;
- mono-Z dilepton and trilepton event selections and lepton-pairing fixes;
- corrected tau flavour construction and robust pileup-weight handling;
- Run-2 mono-Z signal, control-region, pileup, and dataset configurations;
- signal-overlay and combined 2017+2018 plotting utilities;
- HTCondor production helpers with configurable local or EOS output;
- `dnn/` training, scoring, validation, and datacard workflows; and
- `combine/` template, fit-diagnostic, limit, and plotting workflows.

Generated ROOT files, plots, logs, caches, Condor task directories, DNN score
arrays, and Combine fit outputs are intentionally not versioned.  The trained
DNN model checkpoints and compact metadata needed for inference are retained.

Doxygen documentation for C++ code is available [here](http://homepage.iihe.ac.be/~aapopov/hzz2l2nu/doc/). At present it is generated manually by Andrey and might not correspond to the latest version of the code.

## Required packages

Before building this repository, the `MELA` and `MelaAnalytics` packages must have been installed w.r.t. the proper [`LCG`](http://lcginfo.cern.ch) environment.

For `MELA` package
```sh
source env.sh
git clone -b v2.4.2 https://github.com/JHUGen/JHUGenMELA.git
cd JHUGenMELA
git apply ../external/JHUGen_py2to3.patch
cd ..
sh JHUGenMELA/MELA/setup.sh -j 8
cd JHUGenMELA/MELA/data/el9_amd64_gcc12/
chmod +x *.so
cd $HZZ2L2NU_BASE
```

for `MelaAnalytics` package
```sh
git clone -b v2.3 https://github.com/MELALabs/MelaAnalytics
cd MelaAnalytics
source setup.sh
cd $HZZ2L2NU_BASE
```

The versions above intentionally follow the upstream `yujil/hzz2l2nu`
installation.  `JHUGenMELA/` and `MelaAnalytics/` are external checkouts and
are not copied into this repository.

The framework environment is pinned by `env.sh` to `LCG_105` on
`x86_64-el9-gcc12-opt`.  The `pyroothist` dependency is kept as the upstream
Git submodule at `python/hzz/pyroothist` rather than copied into the source
tree.

The DNN utilities additionally use Python 3 with `numpy`, `awkward`, `uproot`,
`matplotlib`, `scikit-learn`, and `torch`.  The statistical workflow in
`combine/` requires a CMSSW environment with the CMS Combine tool installed;
see the [official Combine setup documentation](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/).


## Computing environment and building

At the start of each session, set up the environment with

```sh
. ./env.sh
```

This script also stores the path to the base directory in environment variable `HZZ2L2NU_BASE`, which should then be used in scripts and compiled code to resolve relative paths to auxiliary data files, such as data-driven weights.

Build the package with the following commands:

```sh
git submodule update --init --recursive
cmake -S . -B build
cmake --build build -j "$(nproc)"
```

The warning from CMake about the new version of Boost can be safely ignored. Executable `runHZZanalysis` is put into `$HZZ2L2NU_BASE/bin`, and it is accessible from `$PATH`. To rebuild the package after a change has been introduced to the code, repeat `cmake --build build`. To start the build from scratch, remove the directory `build` and repeat the commands above.

It is also possible to create a program outside of the repository and link it against the shared library of the framework. See [here](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/wikis/shared-library) for documentation.


## Running interactively

Computationally heavy part of the analysis is carried out by program `runHZZanalysis`. Here is an example command to run it interactively:

```sh
runHZZanalysis --config dilepton/2018-ul.yaml \
--ddf /eos/cms/store/group/phys_smp/ZZTo2L2Nu/HZZsample/2018/YAML/DYJetsToLL_PtZ-0To50.yaml \
--analysis DileptonTrees \
--max-events 100000 \
--more-vars \

```
- The first parameter `config` is the path to the master configuration file, such as 2016.yaml. It provides global settings that affect all analyises and all datasets. The path is resolved with the help of FileInPath service. Standard configuration files are located in directory $HZZ2L2NU_BASE/config, which is checked by FileInPath automatically.
- The second parameter `ddf` is the path to a dataset definition file (either a full one or a derived fragment). It provides paths to input files included in the dataset and all dataset-specific configuration parameters.
- The third parameter `analysis` specify which analysis should be executed.
- The fourth parameter `max-events` is the maximal number of events to process.
- The fifth parameter `more-vars`is to add more variables to the processed root.




```sh
runHZZanalysis --analysis <analysis> --help
```

Note that different analyses support different sets of parameters (hence the flag `--analysis` above; without it the help for the default analysis is printed).


## Using batch system

This repository includes scripts to submit jobs to the IIHE batch system and to merge their output.


### Submission

The submission is done with commands like

```sh
prepare_jobs.py --task-dir batch --config 2016.yaml \
  -- $HZZ2L2NU_BASE/config/samples_dilepton_2016.txt -a DileptonTrees
cd batch
big-submission send_jobs.sh
```

The first parameter is the name of the directory in which job scripts and output files will be stored. As with the interactive running, the second parameter is the master configuration file; it is also automatically forwarded to `runHZZanalysis`. The first positional argument is a list of datasets to be processed. Each job will process at maximum one dataset, and large datasets are automatically split among multiple jobs. All subsequent positional arguments are forwarded to the executable without a change. When some of them start with dashes, as in the example above, it is convenient to denote the start of positional arguments with `--` so that they are not interpreted as optional arguments for `prepare_jobs.py`.

By default, `prepare_jobs.py` will run on nodes of the batch system executable `runHZZanalysis`. But it is possible to supply a different executable with option `--prog`.


### Systematic variations

The full list of options supported by `prepare_jobs.py` is available through its help message. One of them is option `--syst` that requests systematic variations to be applied. All supported systematic uncertainties are listed in file [`config/syst.yaml`](config/syst.yaml). For each uncertainty it provides a sequence of masks that defines which datasets are affected by this uncertainty. The masks are checked against the names of the datasets, as specified in [dataset definition files](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions). There are several ways to specify a systematic variation:

* `--syst jec_up` or `--syst jec_down` <br /> These are examples of fully specified variations. They will be propagated to the underlying executable without a change, but only datasets that are affected by those variations (as indicated in `config/syst.yaml`) will be processed.
* `--syst jec` <br /> Jobs for both up and down variations of given type will be created.
* `--syst weights` <br /> Will evaluate all weight-based variations. Only applicable to analyses that produce trees with individual events.
* `--syst all` <br /> Jobs for all possible systematic variations, as well as the nominal configuration, will be created.


### Harvesting

When all jobs have finished (which can be checked with `qstat -u $USER`), their outputs can be merged with

```sh
harvest.py --task-dir batch --config 2016.yaml \
  $HZZ2L2NU_BASE/config/samples_dilepton_2016.txt
```

The will be placed in `batch/merged`. You will probably want to move that directory somewhere and delete the rest of directory `batch`.
## Non-resonant Background Trees

For NRB-tree analysis, there several easy steps produing final trees for templates building. First, take 2016 as an example, prepare the jobs by

```sh
prepare_jobs.py --task-dir nrb2016 --config 2016.yaml -- $HZZ2L2NU_BASE/config/samples_dilepton_2016.txt -a NrbTrees
```

Remember to add `-a NrbTrees` to specify this `NrbTrees` analysis type. You cannot use `--syst` for this analysis since all systematic variations are automatically included. Second, harvest all jobs by

```sh
harvest.py --task-dir nrb2016 --config 2016.yaml $HZZ2L2NU_BASE/config/samples_dilepton_2016.txt
```
Last, produce the final NRB tree for templates building

```sh
nrbTreeHandler -i $HZZ2L2NU_BASE/nrb2016/merged/Data.root -o NRB_weights.root
```
Then, an output root file named `NRB_weights.root` that contains all weights should appear in current directory and you can copy it to the `merged` directory of trees for main analysis. By default, `nrbTreeHandler` will read the `Data.root` in the current directory and compute all weights including systematic variations and add them into different branches. You can also specify the method for the reweighting by

```sh
nrbTreeHandler -i $HZZ2L2NU_BASE/nrb2016/merged/Data.root -a kmethod
```
or
```sh
nrbTreeHandler -i $HZZ2L2NU_BASE/nrb2016/merged/Data.root -a kcorrection
```
or
```sh
nrbTreeHandler -i $HZZ2L2NU_BASE/nrb2016/merged/Data.root -a alphamethod
```
aka. k-method, k-method with ptmiss correction and alpha-method for the estimation. Default is k-method.

## Post-processing with event-based analysis

Instructions below are for the event-based analyses, such as `DileptonTrees`, only. They assume that the harvesting has finished, and the merged files are available in directory `$tree_dir`.

### Plots with comparison between data and simulation

These plots can be produced by running

```sh
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim.yaml --prefix ${tree_dir}/ --output data_sim
```

Event selection and variables to be plotted are described in the configuration file given as the first argument to the script. The [file](config/plot_data_sim.yaml) included in the repository is considered an example, which you would adjust to your needs. Names (or rather name fragments) of input ROOT files are specified in the configuraiton. Each name is _textually_ prepended with the prefix given with the corresponding command line option  (or, alternatively, in the configuration file), which allows to specify the directory containing the files but also a common starting name fragment. In addition to this prefix, multiple standard endings are tried for each file name. The directory in which produced figures will be stored is given by flag `--output`.


## Weights for the photon control region

The following procedure needs to be used to re-compute the weights for the photon CR (in order to get an estimate of the Z+jets contribution in the signal region):

* Run once the `DileptonTrees` analysis (this can be done on data only), using the option `--ptmiss-cut=0`. Harvest the output.
* Run once the `PhotonTrees` analysis (this can be done on data only), changing the config file for NOT applying mean weights (else, it would drop events with a mean weight of 0).
以上两步已经在之前完成。
* Compute nvtx weights:

```sh
compute_instrMET_weights.py -o OUTPUTNAME.root -s nvtx PATH/TO/DILEPTON/TREES PATH/TO/PHOTON/TREES
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_nvtx_2018.root -s nvtx batch_dilepton_2018/merged/ batch_singlephoton_2018/merged/
```

* Move the weights to `data/InstrMetReweighting/`. Erase the previous ones or change the config file.
* Re-run the `PhotonTrees` analysis (the same way as before), making sure that nvtx weights are applied (check that in the config file).
将config/singlephoton/2018-ul.yaml 中的apply_nvtx_reweighting: 改成 True
```sh
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_nvtx_reweight --config singlephoton/${YEAR}-ul-nvtxReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_nvtx_reweight && bash send_jobs.sh); done
6 minutes
harvest.py --task-dir batch_singlephoton_2018_nvtx_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018_data.txt
```
* Compute eta weights, the same way as before:
```sh
compute_instrMET_weights.py -o OUTPUTNAME.root -s eta PATH/TO/DILEPTON/TREES PATH/TO/PHOTON/TREES
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_eta_2018.root -s eta batch_dilepton_2018/merged/ batch_singlephoton_2018_nvtx_reweight/merged/
```

* Re-run the `PhotonTrees` analysis once again, applying the nvtx and eta weights.
将config/singlephoton/2018-ul.yaml 中的apply_eta_reweighting: 也改成 True
```sh
for YEAR in 2018; do prepare_htcondor_jobs.py --task-dir batch_singlephoton_${YEAR}_eta_reweight --config singlephoton/${YEAR}-ul-etaReweight.yaml -- $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_${YEAR}_data.txt --analysis PhotonTrees --more-vars && (cd batch_singlephoton_${YEAR}_eta_reweight && bash send_jobs.sh); done
3 minutes
harvest.py --task-dir batch_singlephoton_2018_eta_reweight/  $HZZ2L2NU_BASE/config/samples_SinglePhoton_InstrMET_2018_data.txt
```
* Compute pT weights. Notice that pT weights also include normalization:

```sh
compute_instrMET_weights.py -o OUTPUTNAME.root -s pt PATH/TO/DILEPTON/TREES PATH/TO/PHOTON/TREES
compute_instrMET_weights.py -o data/InstrMetReweighting/weight_pt_2018.root -s pt batch_dilepton_2018/merged/ batch_singlephoton_2018_eta_reweight/merged/
```

* Compute mass lineshape:

```sh
compute_mass_lineshape.py -o OUTPUTNAME.root -s pt PATH/TO/DILEPTON/TREES
```

* From there, the weights and lineshape can be used in the photon CR analysis (which needs then to be run a fourth time), and mean weights can also be computed.
* To compute mean weights, simply use the script `MeanWeightsComputation.C` in `WeightsAndDatadriven/InstrMET/meanWeights/`. Don't forget to change the names of the input and output files. This is run like a standard ROOT looper:

```sh
root -l
.L MeanWeightsComputation.C++
MeanWeightsComputation t
t.Loop()
.q
```

Don't forget to move the output to the `data/InstrMetReweighting` folder, and to link it correctly from the config file when running the photon trees analysis.


## Producing final (Asimov) templates for the statistical analysis

The following procedure explains the steps needed to have the final (blinded) templates, that are going to be passed to the statistical analysis (to, eventually, make datacards). It assumes that photon reweighting has already been computed, as well as mean weights (if not, see the above section). It assumes you want the Asimov templates, which we need to build by hand using a substraction of data and genuine MET samples in the photon CR as a proxy for the InstrMET pre-fit estimation.

Notice that this procedure needs to be run once with each year if one wants to combine everything in the end.

* Run on dilepton trees, including all systematics (option `syst=all`), using the standard scripts.
* Run on photon trees, including all systematics (option `syst=all`), using the standard scripts (don't forget to allow mean weights; for 2017, don't forget to switch to the correct sim nvtx profiles in the config file).
* Run on NRB trees (see section above), make the final NRB tree and move it with the output of the dilepton trees.
* Prepare to build the templates. At some point, it was possible to do this interactively (using the default options of the script), but with the splitting of the `geq2jets` category, it's not manageable anymore. Therefore, we need to parallelize this and submit jobs on the cluster. These jobs all run the `build_templates.py` script on a given category.
* Prepare the scripts to run on dilepton templates:

```sh
prepare_templates_scripts.py -o templates_dilepton_XXX -s dilepton -y YEAR PATH/TO/DILEPTON/TREES
```

* This should have created a bunch of config files inside a directory `templates/scripts/templates_dilepton_XXX`. Go to that directory and launch the jobs (there should be 9 jobs):

```sh
cd templates/scripts/templates_dilepton_XXX
big-submission submit.sh
```

* Once the jobs have finished, you should have 9 separate output files inside the `templates` folder. Just combine them using hadd. CAUTION: this, for whatever reason, does not work on the current HZZ2l2nu fwk infrastructure; I found a work-around by doing a `cmsenv` from a recent CMSSW release (like `10_2_22`).

```sh
cd ${HZZ2L2NU_BASE}/templates
hadd templates_dilepton_XXX.root templates_dilepton_XXX_*.root
rm templates_dilepton_XXX_*.root
```

* Do NOT look at the data in these templates: these are filtered with final analysis cuts, and we don't want to look at them as we are still unblinded.
* Do the same for photon templates (or run it directly interactively).
* Build an Asimov template that we will use instead of the template for the signal region:

```sh
build_asimov_templates.py -w ../data/InstrMetReweighting/meanWeights_YEAR.root -o templates_dilepton_XXX_Asimov.root templates_dilepton_XXX.root templates_photon_XXX.root
```

* From there, we want to prepare the templates for the statistical analysis, by splitting the CR into one CR per mT bin, each of these bins containing one process. This is done with the script `WeightsAndDatadriven/makePhotonTemplates.C`: make sure that you change the names in that script (inputs and outputs: lines 2, 3 and 14). (Make sure you are using Asimov templates for data and not real ones!)

```sh
cd ${HZZ2L2NU_BASE}/WeightsAndDatadriven
root -l -q 'makePhotonTemplates.C(false)'
root -l -q 'makePhotonTemplates.C(true)'
```

* Once again, go to a build that allows to `hadd` properly (like you CMSSW release) for the next command:

```sh
hadd -f templates_dilepton_XXX_Asimov_withInstrMET.root temporary_SR_bin_*.root
rm temporary_SR_bin*.root
```

* The rest is handled separately in the repository for the statistical analysis. The templates that we will use there are the many different templates for the CR (one per bin) and the template for CR that we just built.



## Interacting with the IIHE cluster

1. Check my jobs

   ```sh
   qstat -u $USER
   ```

2. Delete jobs

   ```sh
   for j in $(qselect -u $USER); do timeout 3 qdel -a $j; done
   ```

3. Find jobs name and info

   ```sh
   qstat -f $(qselect -u $USER -s EHQRTW) | grep "Job Id\|Job_Name\|resources_used.walltime"
   ```

4. Go interactively on the express queue for some tests

   ```sh
   qsub -I -q express
   ```

5. Peek at your running jobs

   - Connect to an `m*` machine with option `-A`.
   - Do `qpeek <jobID>`.

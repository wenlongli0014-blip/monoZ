# DNN workflow

This directory contains the reusable source, configuration, trained checkpoints,
and datacards for the mono-Z DNN studies.  Generated ROOT files, plots, logs,
intermediate arrays, and unpublished numerical summaries are intentionally not
versioned.

The main entry points are:

- `train_dnn_score.py` for training;
- `add_dnn_score_to_root.py` and `make_dnn_shapes.py` for inference and shapes;
- `make_dnn_datacard.py` for datacard construction;
- the `run_*.sh` scripts for the year-specific and combined workflows.

The Python environment needs PyTorch, NumPy, Awkward Array, uproot,
scikit-learn, Matplotlib, and PyROOT.  Statistical inference additionally needs
the CMS Combine environment described in `../README.md`.

Before running a workflow, review the input paths in the corresponding JSON
configuration under `2017/`, `2018/`, or `combined_2017_2018/` and adapt them to
the local storage layout.

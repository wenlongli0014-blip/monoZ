# Statistical-analysis workflow

This directory keeps the reusable scripts and hand-written inputs for building
mono-Z shapes, datacards, Asimov datasets, fits, and limit scans with CMS
Combine.  Generated workspaces, ROOT files, plots, logs, tables, and unpublished
thesis text are intentionally kept outside version control.

The main groups of scripts are:

- `make_shapes_*.py` and `make_pseudo_asimov_data.py` for statistical inputs;
- `make_datacards_*.py` for datacard construction;
- `add_manual_shape_uncertainties.py` for the manual shape variations;
- `make_yield_tables.py` and the plotting scripts in the subdirectories for
  local validation and presentation.

Run these tools inside a compatible CMS Combine environment.  The recommended
Combine installation link and the common analysis dependencies are documented
in `../README.md`.  Input paths and analysis choices are site-specific and
should be reviewed before executing a workflow.

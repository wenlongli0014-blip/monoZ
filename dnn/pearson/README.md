# HZZ DNN Pearson correlation audit

This directory contains train-only Pearson correlation diagnostics for the
current 2017, 2018, and combined 2017+2018 DNN inputs.

See `REPORT.md` for the main numerical findings from the generated matrices.

Run the complete audit with:

```bash
python3 dnn/pearson/compute_pearson.py
```

The script deliberately follows the current DNN training definitions:

- explicit candidate features from each DNN JSON config;
- SR selection `lepton_cat != 2 && ptmiss >= 100`;
- a reproducible 50% class-stratified training split with seed `12345`;
- non-negative event weights clipped at the configured value (`100`);
- pairwise removal of non-finite values and the `-9999` missing sentinel.

For each dataset, outputs are split into `candidate_features` (all 36 inputs)
and `selected_top10` (the inputs used by the saved model). Only physically
weighted Pearson results are retained. Each directory has:

- `pearson_{all,signal,background}_weighted.png`: annotated heatmaps;
- matching `.csv` files: exact correlation matrices;
- `*_effective_events.csv`: pairwise effective weighted event counts;
- `*_pair_event_count.csv`: pairwise raw valid event counts;
- `pearson_pairs.csv`: every feature pair, sorted by its largest absolute
  correlation across the six matrices;
- `strong_pairs.csv`: the subset with maximum absolute correlation at least
  `0.80`;
- `pearson_matrices.npz`: all matrices in machine-readable NumPy format.

Each dataset also has `feature_provenance.csv`, which classifies inputs as a
direct ROOT scalar, an indexed element extracted from a ROOT vector branch, or
a feature derived during DNN input construction.

The `all` matrix is useful for comparison but should not be used alone for
feature removal: class mixing can create or hide correlations. Signal and
background matrices are therefore the primary diagnostics.

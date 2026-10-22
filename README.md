# IsoCred

This repository has the code and processed data used in our comparison of CAGE-supported promoters in HepG2 and WTC11. You can check the reported numbers and rerun predictions from the saved models without downloading the raw sequencing data.

In HepG2, the tested hierarchical MAP model did not improve on the strongest supervised baseline. We then compared models with matched inputs and examined detection opportunity, using development data after the original holdout had been evaluated. In WTC11, the separately specified common10 models reversed the pooled RF–elastic-net ordering. The files also retain the missing predictions from a failed elastic-net fold.

## Quick check

Use Python 3.11 (the analysis used 3.11.16). Run from this directory:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py verify
python reproduce.py replay
```

`verify` checks file hashes and recomputes stored point estimates from group-level labels and predictions. `replay` also runs the included fixed checkpoints without fitting: HepG2 final MAP/RF and baselines, post-holdout RF14/EN14/logistic controls, and WTC11 RF10/EN10. Both only write a new `replay-output/` directory. Read the output JSON for the exact checks and limitations. Only load the checkpoints shipped with this release: joblib files can execute code when loaded.

## What is here

The files include group features and labels, dependency and fold identifiers, eligibility, predictions, tuning settings, saved models and result tables. The model, baseline and endpoint code used in the study is included too. `DATA.md` explains the inputs; `provenance/source_to_release.json` records where each copied file came from and any edits needed to make paths portable.

A few places to start:

- For the original MAP comparison and its mechanism controls, see the [development metrics](work/ISO-001A/output/results/ISO-001A_metrics.tsv).
- For the original HepG2 holdout, see the [main metrics](output/ISO-002B/final_metrics/evaluation/metrics.tsv) and [paired intervals](output/ISO-002B/final_metrics/evaluation/cluster_intervals.tsv).
- For matched inputs and detection opportunity after the holdout, see the [development contrasts](work/ISO-002C/evaluation/results_A/contrast_CI.tsv) and [common-support comparisons](work/ISO-002C/evaluation/results_A/paired_common_support_models.tsv).
- For WTC11, the [five-method prediction table](output/ISO-002D/predictions/WTC11_all_five_methods.parquet) can be reused directly. The [complete results table (gzip)](work/ISO-002D/evaluation_adapter/run/external/scientific_numbers.tsv.gz) includes the primary comparison and every sensitivity endpoint.

The quick checks start from processed group tables. Scripts under `work/` document the earlier preprocessing, fitting and 2,000-draw component bootstrap. Those scripts retain their original protocol and hash checks, so they need the documented inputs rather than a single raw-data command. Discovery, training and the full bootstrap were not rerun for this release.

For raw reads and CTSS files, use the source accessions in `DATA.md`. Raw data, third-party software distributions, the manuscript and artwork are not bundled here.

The measured outcome is promoter support. It does not establish the function or full structure of a transcript. WTC11 tests the common10 models under changes in cell system, candidate generation and alignment; it does not test the original RF71 model or require unseen genes. Davies/A375 and other reserved external outcomes are not part of this release.

## License

See `LICENSE` and `NOTICE` for the code/data distinction and upstream-resource terms.

# IsoCred

Code and processed data for our study of CAGE-supported promoters in HepG2 and WTC11. The commands below let you check the reported results and rerun predictions from saved models without downloading the raw sequencing data.

The tested hierarchical MAP model did not improve on the strongest supervised baseline in HepG2. We followed up with matched-input and detection-opportunity comparisons on development data, after the original holdout had been evaluated. In WTC11, the separately specified common10 models reversed the pooled RF–elastic-net ordering. Predictions from the failed elastic-net fold remain missing.

## Quick check

Use Python 3.11 on Linux; the analysis used 3.11.16. Run from this directory:

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python reproduce.py verify
python reproduce.py replay
```

`verify` checks file hashes and recalculates point estimates from the included labels and predictions. `replay` reruns saved-model predictions for the HepG2 final comparison, the later RF14/EN14/logistic controls, and WTC11 RF10/EN10. It does not fit models. Results go to `replay-output/verify/` and `replay-output/replay/`; existing results are never overwritten. Use `--output <new-directory>` to run a check again.

Load the checkpoints only if you trust this release; joblib deserialization can execute code.

## What is here

The files include group features and labels, fold and dependency identifiers, predictions, tuning settings, saved models, result tables, and the analysis code. [DATA.md](DATA.md) explains the inputs. `provenance/` records source files, software versions and checksums.

A few places to start:

- For the original MAP comparison and its mechanism controls, see the [development metrics](work/ISO-001A/output/results/ISO-001A_metrics.tsv).
- For the original HepG2 holdout, see the [main metrics](output/ISO-002B/final_metrics/evaluation/metrics.tsv) and [paired intervals](output/ISO-002B/final_metrics/evaluation/cluster_intervals.tsv).
- For matched inputs and detection opportunity after the holdout, see the [development contrasts](work/ISO-002C/evaluation/results_A/contrast_CI.tsv) and [common-support comparisons](work/ISO-002C/evaluation/results_A/paired_common_support_models.tsv).
- For WTC11, the [five-method prediction table](output/ISO-002D/predictions/WTC11_all_five_methods.parquet) can be reused directly. The [complete results table (gzip)](work/ISO-002D/evaluation_adapter/run/external/scientific_numbers.tsv.gz) includes the primary comparison and every sensitivity endpoint.

The quick checks start from processed group tables. Scripts under `work/` record the earlier preprocessing, fitting and 2,000-draw component bootstrap. They retain the original protocol and hash checks and need the corresponding inputs; they are not a one-command pipeline from raw data. Discovery, training and the full bootstrap were not rerun for this release.

The plans in `contracts/` are historical records. Their dates, management notes and earlier submission plans describe that stage of the study. Scripts and model records identify the original protocol hashes; the source inventory records the copies included here and any path edits.

For raw reads and CTSS files, use the source accessions in `DATA.md`.

The measured outcome is promoter support, not the function or full structure of a transcript. WTC11 tests the common10 models under changes in cell system, candidate generation and alignment. It is not an external test of RF71, and the two datasets can contain the same genes.

## License

See `LICENSE` and `NOTICE` for the code/data distinction and upstream-resource terms.

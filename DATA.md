# Data and scope

All tables are keyed by `promoter_group_id`. Grouping preserves a maximum 50-bp all-member TSS envelope. A group supplies one target and prediction; it is not an independently validated full transcript. Exact boundary ties are retained for top-fraction enrichment, so actual selected coverage can exceed the nominal fraction. NaN means unavailable/unscored, never an observed zero.

## Public inputs

* HepG2 long-read inputs: ENCODE ENCFF427JDY and ENCFF589SMB. The historical calls used TALON 6.0.1, FLAIR 3.0.1, IsoQuant 4.0.0, SQANTI3 6.0.2, minimap2 2.31 and SAMtools 1.24, with GRCh38/GENCODE v29.
* HepG2 CAGE: FANTOM5 hg38 CTSS, CNhs12328/DRZ012937, CNhs12329/DRZ012938 and CNhs12330/DRZ012939. Positive requires raw signal >0 within inclusive same-strand ±50 bp of any exact member start in at least two of three biological replicates. S2 additionally requires positions within 5 bp. The group labels already unlocked for these analyses are included; raw CTSS is not.
* WTC11 long-read neutral measurements: ENCODE ENCFF166IXJ (source ENCFF563QZR; bio1 ENCBS944CBA) and ENCFF708JUG (source ENCFF245IPA; bio3 ENCBS474NOC). Original archive models were `MortazaviLab/Challenge1/iso_detect_ref_talon_lapa_long/human_WTC11_cDNA_PacBio/models.gtf.gz` from the [LRGASP archive](https://cgl.gi.ucsc.edu/data/LRGASP/submissions/MortazaviLab.tar.gz) (archive SHA-256 `15b8cfcd3820fab01f6fe5af3cbd8841b55a834fead8d6bad6a95d8b7f6bb0f7`); original model coordinates, not corrected replacements, define the candidates. The archive declared TALON 4.0.0 and LAPA 0.0.2; the external alignment used minimap2 2.15 splice/-C5 rather than the internal minimap2 2.31 splice:hq setting. Exact archival execution commit/config remains unbound. The archive source identity and raw-to-group reconstruction are outside this processed-table replay; do not substitute a newly generated transcript set and call it exact reproduction.
* WTC11 ssCAGE: GEO GSM5625727, GSM5625728, GSM5625729 (bio1 technical libraries pooled by coordinate) and GSM5625730 (bio3). Primary support requires both biological pools (2/2). This is two biological pools, one donor, not four independent biological replicates. All three technical-pair sensitivities and the 5-bp endpoint remain in the frozen results. GSM5625731/32 were excluded under the fixed sample contract.

Retrieve ENCODE files at `https://www.encodeproject.org/files/<accession>/`; GEO series/sample pages at `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>`; FANTOM5 resources at `https://fantom.gsc.riken.jp/5/`. Upstream resources retain their own licenses and access conditions.

## Analysis boundaries

* Original HepG2 development: 30,487 eligible groups, five outer folds with dependency-component purging and two fixed inner hash partitions. EN71 fold 3 failed and remains missing (24,223 scored groups).
* Original HepG2 holdout: 11,602 assigned groups, 10,464 eligible, 773 dependency components. Model predictions were fixed before outcome exposure. MAP-minus-RF top-decile enrichment is negative; the original algorithm STOP is preserved.
* Post-holdout matched-input/detection-opportunity analysis: development-only, not another untouched HepG2 test. RF14/EN14 and three logistic controls are included with their original failures, parameters and common-support analyses.
* Common10: ten named features; trained only in eligible HepG2 development groups. WTC11 has 35,560 structurally measured groups and 34,150 primary-eligible groups in 3,997 components. RF–EN ranking reverses for top-decile enrichment; the AUC interval spanning zero does not establish equivalence.

Data folders retain historical relative names to keep the source crosswalk unambiguous. The four HepG2 group feature tables hold represented candidate-core columns, evidence summaries, baselines and dependency/split metadata. common10 excludes read-count/shared-read descriptors unavailable under the fixed external representation. Training-only median imputation, missing indicators and robust scaling are stored in each checkpoint. Full-tie and fractional-boundary sensitivity results, subgroup coverage and uncertainty are kept in the result tables. Original 2,000-draw paired component-bootstrap code is supplied but not rerun by the quick check.

Excluded: Davies/A375 outcomes, globally reserved external outcomes, exposure exclusions, raw sequencing reads, raw CTSS, third-party implementation copies, machine logs, credentials, author details, manuscripts and images.

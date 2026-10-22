# Data and scope

Group-level tables use `promoter_group_id` as the key. Within each group, the earliest and latest member TSS are at most 50 bp apart. Each group has one target and one prediction. This measures promoter support rather than verifying every transcript structure. All exact ties at a selection boundary are retained, so selected coverage can exceed the nominal fraction. NaN means unavailable or unscored; it is not an observed zero.

## Public inputs

* HepG2 long-read inputs: ENCODE ENCFF427JDY and ENCFF589SMB. The historical calls used TALON 6.0.1, FLAIR 3.0.1, IsoQuant 4.0.0, SQANTI3 6.0.2, minimap2 2.31 and SAMtools 1.24, with GRCh38/GENCODE v29.
* HepG2 CAGE: FANTOM5 hg38 CTSS, CNhs12328/DRZ012937, CNhs12329/DRZ012938 and CNhs12330/DRZ012939. A group is positive when raw signal >0 occurs within inclusive same-strand ±50 bp of any exact member start in at least two of three biological replicates. S2 additionally requires the supporting positions to be within 5 bp. The processed group labels are included; original CTSS files remain with the provider.
* WTC11 long-read measurements: ENCODE ENCFF166IXJ (source ENCFF563QZR; bio1 ENCBS944CBA) and ENCFF708JUG (source ENCFF245IPA; bio3 ENCBS474NOC). Candidate coordinates come from the original `MortazaviLab/Challenge1/iso_detect_ref_talon_lapa_long/human_WTC11_cDNA_PacBio/models.gtf.gz` in the [LRGASP archive](https://cgl.gi.ucsc.edu/data/LRGASP/submissions/MortazaviLab.tar.gz), not corrected replacements. The archive SHA-256 is `15b8cfcd3820fab01f6fe5af3cbd8841b55a834fead8d6bad6a95d8b7f6bb0f7`. It declared TALON 4.0.0 and LAPA 0.0.2. WTC11 alignment used minimap2 2.15 splice/-C5; HepG2 used minimap2 2.31 splice:hq. The exact generator commit and configuration for the archived run could not be established. The checks here start from processed tables rather than rebuilding groups from raw data. A newly generated transcript set would not be an exact replacement for the archived candidates.
* WTC11 ssCAGE: GEO GSM5625727, GSM5625728, GSM5625729 (bio1 technical libraries pooled by coordinate) and GSM5625730 (bio3). Primary support requires both biological pools (2/2). This is two biological pools, one donor, not four independent biological replicates. All three technical-pair sensitivities and the 5-bp endpoint remain in the frozen results. GSM5625731/32 were excluded under the fixed sample contract.

Retrieve ENCODE files at `https://www.encodeproject.org/files/<accession>/`; GEO series/sample pages at `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>`; FANTOM5 resources at `https://fantom.gsc.riken.jp/5/`. Upstream resources retain their own licenses and access conditions.

## Analysis boundaries

* Original HepG2 development: 30,487 eligible groups, five outer folds with dependency-component purging and two fixed inner hash partitions. EN71 fold 3 failed and remains missing (24,223 scored groups).
* Original HepG2 holdout: 11,602 assigned groups, 10,464 eligible, 773 dependency components. Model predictions were fixed before outcomes were read. MAP-minus-RF top-decile enrichment is negative. Development of the original MAP method stopped after it did not meet the prespecified comparison criteria.
* Post-holdout matched-input/detection-opportunity analysis: development-only, not another untouched HepG2 test. RF14/EN14 and three logistic controls are included with their original failures, parameters and common-support analyses.
* Common10: ten named features; trained only in eligible HepG2 development groups. WTC11 has 35,560 structurally measured groups and 34,150 primary-eligible groups in 3,997 components. RF–EN ranking reverses for top-decile enrichment; the AUC interval spanning zero does not establish equivalence.

Folders retain their original relative names so that file references still match the analysis scripts. The four HepG2 group feature tables contain the represented core features, evidence summaries, baselines and dependency/split metadata. common10 excludes read-count and shared-read descriptors that could not be reconstructed with the same meaning in WTC11. Each checkpoint includes the imputation, missing indicators and robust scaling fitted on its training data. Result tables retain full-tie and fractional-boundary sensitivities, subgroup coverage and uncertainty. The 2,000-draw paired component-bootstrap code is included, but the quick check does not rerun it.

Large upstream sequencing and CTSS files are linked rather than bundled. See `NOTICE` for attribution and upstream terms.

#!/usr/bin/env python3
"""Frozen Phase-1 development OOF statistics. Never fits or reads final/external outcomes."""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import argparse
import hashlib
import itertools
import json
import math
import platform
import resource
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from threadpoolctl import threadpool_limits
threadpool_limits(1)
pa.set_cpu_count(1)
pa.set_io_thread_count(1)

ROOT = Path(".")
HERE = ROOT / "work/ISO-002C/evaluation"
CONTRACT = "contracts/ISO-002C_phase1_comparison_contract.yaml"
CONTRACT_SHA = "b1c0099629f6fb67d9b620951fabdd6994aead597b93c520d9ffb900eae5829d"
SEED, DRAWS, CPU_LIMIT = 20260903, 2000, 3600
EVIDENCE = "POST_HOC_SECONDARY_AFTER_FINAL_EXPOSURE_DEVELOPMENT_GROUPED_OOF"
OLD = "work/ISO-001A/output/results/ISO-001A_development_OOF.parquet"
META = "work/ISO-001A/output/data/group_metadata_label_blind.parquet"
LABELS = "work/ISO-001A/output/data/development_E4_v1.parquet"
OLD_MANIFEST = "logs/ISO-002B_checkpoint_manifest.json"
OLD_MANIFEST_SHA = "1e02ec1e89e0f8dcaaed623707a3e666b9a8ac74011bd1190428d2e384cdf511"
FIXED_HASHES = {
    CONTRACT: CONTRACT_SHA,
    OLD_MANIFEST: OLD_MANIFEST_SHA,
    OLD: "661b9089ce435cb33fc7413822860bab9b7048e3cb76f68d8e2d05bf83c6e994",
    META: "2df9945cf0f330b17de9787c72881345d5ede921711f64884e2314269515c346",
    LABELS: "ef9bd0b3fcb6a9fca6032d06849162e696aa0ad2651fe6d53113783e4de958f0",
}
NEW = ("RF14", "EN14", "L2_READ", "L2_READ_ANNOTATION", "L2_READ_ANNOTATION_OPPORTUNITY")
ALIASES = {
    "RF71": "RANDOM_FOREST", "EN71": "ELASTIC_NET", "MAP14": "ISOCRED_MAP",
    "READ_UNION": "READ_UNION", "B0": "B0", "B2": "B2",
    "SOURCE_AWARE": "SOURCE_AWARE", "PSQAN_STYLE": "PSQAN_STYLE",
    "SQANTI3_RULE": "SQANTI3_RULE",
}
MODELS = tuple(ALIASES) + NEW
PROBABILITIES = {"RF71", "EN71", "MAP14", *NEW}
CONTRASTS = (
    ("RF14", "MAP14"), ("EN14", "MAP14"), ("RF71", "RF14"),
    ("EN71", "EN14"), ("L2_READ", "READ_UNION"),
    ("L2_READ_ANNOTATION", "L2_READ"),
    ("L2_READ_ANNOTATION_OPPORTUNITY", "L2_READ_ANNOTATION"),
    ("RF71", "L2_READ_ANNOTATION_OPPORTUNITY"),
)
LEVELS = {
    "read": ("0", "1-4", "5-19", "20-99", ">=100"),
    "annotation": ("novel-only", "mixed", "known-only"),
    "sharing": ("1", "2-4", ">=5"),
    "window": ("101", "102-125", ">125"),
    "monoexon": ("none", "mixed", "all"),
    "caller_TALON": ("absent", "partial", "all"),
    "caller_FLAIR": ("absent", "partial", "all"),
    "caller_IsoQuant": ("absent", "partial", "all"),
}
JOINT = tuple("|".join(c) for c in itertools.product(
    LEVELS["read"], LEVELS["annotation"], LEVELS["sharing"]))
Q = (0.05, 0.10, 0.20)
MODES = ("full_ties", "fractional_boundary")
LOG_EPS = 1e-12

def require(test, message):
    if not bool(test):
        raise ValueError(message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def clean_json(x):
    if isinstance(x, dict):
        return {str(k): clean_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in x]
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(x) else None
    return x

def dump_new(path, obj):
    with Path(path).open("x") as f:
        json.dump(clean_json(obj), f, indent=2, allow_nan=False)
        f.write("\n")

def file_record(rel, role):
    p = ROOT / rel
    require(p.is_file() and not p.is_symlink(), "Missing file or symlink: " + rel)
    require(p.resolve().is_relative_to(ROOT.resolve()), "Input path escape: " + rel)
    return {"path": rel, "role": role, "bytes": p.stat().st_size, "sha256": sha(p)}

def new_paths():
    return [f"output/ISO-002C/results/{m}_OOF.parquet" for m in NEW]

def training_receipts():
    return [f"output/ISO-002C/receipts/{m}_training_complete.json" for m in NEW]

def whitelist():
    return tuple(FIXED_HASHES) + tuple(new_paths()) + tuple(training_receipts())

def check_code_freeze():
    p = HERE / "receipts/code_freeze.json"
    f = json.loads(p.read_text())
    require(f["freeze_status"] == "frozen", "Evaluation code not frozen")
    required = {
        "work/ISO-002C/evaluation/evaluate_phase1.py",
        "work/ISO-002C/evaluation/tests/test_evaluate_phase1.py",
        "work/ISO-002C/evaluation/README.md",
    }
    entries = {r["path"]: r for r in f["files"]}
    require(required.issubset(entries), "Code freeze lacks evaluator/tests/interface")
    for rel, r in entries.items():
        observed = file_record(rel, r["role"])
        require(observed["bytes"] == r["bytes"] and observed["sha256"] == r["sha256"],
                "Code freeze drift: " + rel)
    return sha(p)

def bind_inputs():
    """Hash-only binding: no Parquet outcome/score deserialization or statistics."""
    code_hash = check_code_freeze()
    dest = HERE / "receipts/input_binding.json"
    require(not dest.exists(), "Refusing to replace existing input binding")
    records = []
    for rel in whitelist():
        r = file_record(rel, "DEVELOPMENT_INPUT" if rel.endswith(".parquet") else "GOVERNANCE")
        if rel in FIXED_HASHES:
            require(r["sha256"] == FIXED_HASHES[rel], "Historical/contract input drift: " + rel)
        records.append(r)
    for rel in training_receipts():
        x = json.loads((ROOT / rel).read_text())
        require(x["contract_sha256"] == CONTRACT_SHA, "Training receipt contract drift")
        require(x["status"] in ("COMPLETE", "COMPLETE_WITH_FAILURES_RETAINED"),
                "Model training not terminal")
        require(x["eligible"] == 30487, "Training eligibility changed")
    dump_new(dest, {
        "schema_version": 1, "status": "FROZEN_INPUT_HASH_BINDING_BEFORE_REAL_EVALUATION",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "code_freeze_sha256": code_hash, "contract_sha256": CONTRACT_SHA,
        "files": records, "parquet_content_read_at_binding": False,
        "no_final_external_or_raw_E4_inputs": True,
    })
    print(json.dumps({"status": "INPUTS_HASH_BOUND_NOT_EVALUATED", "sha256": sha(dest)}))

def verify_bound_inputs():
    code_hash = check_code_freeze()
    p = HERE / "receipts/input_binding.json"
    b = json.loads(p.read_text())
    require(b["code_freeze_sha256"] == code_hash, "Input binding code freeze mismatch")
    require(b["contract_sha256"] == CONTRACT_SHA, "Bound contract mismatch")
    entries = {r["path"]: r for r in b["files"]}
    require(len(entries) == len(b["files"]) and set(entries) == set(whitelist()),
            "Input registry not equal to exact whitelist")
    for rel, r in entries.items():
        o = file_record(rel, r["role"])
        require(o["sha256"] == r["sha256"] and o["bytes"] == r["bytes"], "Input drift: " + rel)
        if rel in FIXED_HASHES:
            require(o["sha256"] == FIXED_HASHES[rel], "Fixed input drift")
    c = json.loads((ROOT / CONTRACT).read_text())
    require(c["primary_supplementary_contrasts"] == [list(v) for v in CONTRASTS],
            "Declared contrast family changed")
    require(c["statistics"]["bootstrap_draws"] == DRAWS and
            c["statistics"]["bootstrap_seed"] == SEED, "Bootstrap contract changed")
    return b, sha(p)

def indexed(x):
    require("promoter_group_id" in x, "Missing group key")
    require(x.promoter_group_id.notna().all() and not x.promoter_group_id.duplicated().any(),
            "Null/duplicate group keys")
    return x.set_index("promoter_group_id")

def merge_scores(old, additions):
    """ID-based exact-cohort join; never align by row order or impute missing scores."""
    d = indexed(old).copy()
    for model, incoming in additions.items():
        require(set(incoming.columns) == {"promoter_group_id", model}, "Unexpected new OOF schema")
        z = indexed(incoming)
        require(set(z.index) == set(d.index), "Different eligible cohort for " + model)
        d[model] = z.loc[d.index, model]
    return d

def load_data():
    d = merge_scores(pd.read_parquet(ROOT / OLD),
                     {m: pd.read_parquet(ROOT / p) for m, p in zip(NEW, new_paths())})
    require(len(d) == 30487, "Eligible universe must remain 30487")
    d = d.rename(columns={v: k for k, v in ALIASES.items() if v != k})
    meta = indexed(pd.read_parquet(ROOT / META)).loc[d.index]
    lab = indexed(pd.read_parquet(ROOT / LABELS))
    require(set(lab.index[lab.evaluation_eligible]) == set(d.index), "Label eligible universe drift")
    lab = lab.loc[d.index]
    require(lab.evaluation_eligible.all() and
            np.array_equal(lab.primary_target.to_numpy(), d.target.to_numpy()),
            "Historical OOF target differs from frozen development labels")
    require(meta.split.eq("DEVELOPMENT_POOL").all() and not meta.touches_final_dependency.any(),
            "Non-development or final-touching group")
    for col in ("dependency_component_id", "inner_fold"):
        require(np.array_equal(d[col].to_numpy(), meta[col].to_numpy()), "Metadata mismatch " + col)
    require(set(d.inner_fold) == set(range(5)), "Frozen folds missing")
    require(set(d.target) == {0, 1}, "Positive target must be binary class 1")
    require(d.dependency_component_id.notna().all(), "Missing dependence component")
    require(d.EN71.isna().sum() == 6264 and
            np.array_equal(d.EN71.isna().to_numpy(), d.inner_fold.eq(3).to_numpy()),
            "Historical EN71 missing-fold pattern changed")
    for m in MODELS:
        score = d[m].to_numpy(dtype=float)
        require(not np.isinf(score).any(), "Infinite scores " + m)
        if m in PROBABILITIES:
            require(((score[np.isfinite(score)] >= 0) & (score[np.isfinite(score)] <= 1)).all(),
                    "Invalid positive-class probabilities " + m)
    for m, p in zip(NEW, training_receipts()):
        r = json.loads((ROOT / p).read_text())
        require(int(d[m].notna().sum()) == r["scored"], "Receipt scoring count mismatch " + m)
        for _, z in d.groupby("inner_fold"):
            require(z[m].notna().all() or z[m].isna().all(), "Partial scoring within outer fold " + m)
    return d, meta

class Ranked:
    """Tie-aggregated weighted rank statistics; no group ID participates in ranking."""
    def __init__(self, y, score, probability=False):
        y, score = np.asarray(y, dtype=float), np.asarray(score, dtype=float)
        require(y.shape == score.shape and y.ndim == 1, "Rank input dimensions")
        require(np.isin(y, [0, 1]).all(), "Rank target is not binary 0/1")
        require(np.isfinite(score).all(), "Ranked input must use explicitly finite common support")
        if probability:
            require(((score >= 0) & (score <= 1)).all(), "Probability outside [0,1]")
        self.order = np.argsort(-score, kind="stable")
        self.y, self.s = y[self.order], score[self.order]
        self.start = np.r_[0, np.flatnonzero(self.s[:-1] != self.s[1:]) + 1] if len(y) else np.array([], dtype=int)
        self.probability = probability
        self.error = (self.s - self.y) ** 2
        p = np.clip(self.s, LOG_EPS, 1 - LOG_EPS)
        self.loss = -(self.y * np.log(p) + (1 - self.y) * np.log1p(-p))

    def calculate(self, weights=None):
        w = np.ones(len(self.y)) if weights is None else np.asarray(weights, dtype=float)[self.order]
        require(len(w) == len(self.y) and np.isfinite(w).all() and (w >= 0).all(), "Invalid bootstrap weights")
        n, pos = float(w.sum()), float(w @ self.y)
        out = {"weighted_n": n, "weighted_positives": pos,
               "prevalence": pos / n if n else np.nan,
               "prevalence_ceiling": n / pos if pos else np.nan,
               "AUC": np.nan, "within_pairs": pos * (n - pos),
               "Brier": np.nan, "log_loss": np.nan}
        if not n:
            total = positive = cs = cp = np.array([])
        else:
            total = np.add.reduceat(w, self.start)
            positive = np.add.reduceat(w * self.y, self.start)
            cs, cp = np.cumsum(total), np.cumsum(positive)
            negative = total - positive
            if 0 < pos < n:
                out["AUC"] = float((positive @ (n - pos - np.cumsum(negative) + .5 * negative)) / (pos * (n - pos)))
            if self.probability:
                out["Brier"], out["log_loss"] = float(w @ self.error / n), float(w @ self.loss / n)
        for q in Q:
            for mode in MODES:
                key = f"{mode}_{int(q * 100):02d}"
                if not n:
                    vals = (np.nan,) * 5
                else:
                    k = min(int(np.searchsorted(cs, q * n, side="left")), len(cs) - 1)
                    before, before_pos = (cs[k-1], cp[k-1]) if k else (0., 0.)
                    factor = 1. if mode == "full_ties" else (q * n - before) / total[k]
                    selected, selected_pos = before + factor * total[k], before_pos + factor * positive[k]
                    ppv = selected_pos / selected
                    vals = (ppv / (pos/n) if pos else np.nan, ppv, selected/n, selected,
                            float(self.s[self.start[k]]))
                for metric, value in zip(("enrichment", "PPV", "selection_coverage", "selected_weight", "cutoff"), vals):
                    out[f"{key}_{metric}"] = value
        return out

class Prepared:
    def __init__(self, y, score, folds, probability=False):
        self.rank = Ranked(y, score, probability)
        self.fold_indices = [np.flatnonzero(np.asarray(folds) == f) for f in range(5)]
        self.fold_ranks = [Ranked(np.asarray(y)[ix], np.asarray(score)[ix], probability) for ix in self.fold_indices]
    def calculate(self, w=None, folds=False):
        w = np.ones(len(self.rank.y)) if w is None else np.asarray(w)
        out = self.rank.calculate(w)
        ff = [r.calculate(w[ix]) for ix, r in zip(self.fold_indices, self.fold_ranks)]
        den = sum(v["within_pairs"] for v in ff)
        out["within_fold_pair_weighted_AUC"] = (
            sum(v["AUC"] * v["within_pairs"] for v in ff if v["within_pairs"] > 0) / den
            if den else np.nan)
        return (out, ff) if folds else out

def make_strata(meta):
    r = meta.read_support_union.to_numpy()
    c = meta.candidate_count.to_numpy()
    win = meta.window_opportunity_bp.to_numpy()
    require((r >= 0).all() and (r == np.floor(r)).all(), "Invalid read count")
    require((c >= 1).all() and (c == np.floor(c)).all(), "Invalid sharing count")
    require((win >= 101).all() and (win == np.floor(win)).all(), "Invalid window opportunity")
    out = {
        "read": np.select([r == 0, r <= 4, r <= 19, r <= 99], LEVELS["read"][:-1], default=">=100"),
        "sharing": np.select([c == 1, c <= 4], LEVELS["sharing"][:-1], default=">=5"),
        "window": np.select([win == 101, win <= 125], LEVELS["window"][:-1], default=">125"),
    }
    for axis, col, zero, one, other in (
        ("annotation", "known_fraction", "novel-only", "known-only", "mixed"),
        ("monoexon", "monoexon_fraction", "none", "all", "mixed"),
        *[(f"caller_{c}", f"caller_{c}_candidate_fraction", "absent", "all", "partial")
          for c in ("TALON", "FLAIR", "IsoQuant")],
    ):
        v = meta[col].to_numpy()
        require(np.isfinite(v).all() and ((v >= 0) & (v <= 1)).all(), "Invalid strata fraction " + col)
        out[axis] = np.select([v == 0, v == 1], [zero, one], default=other)
    out["joint_read_annotation_sharing"] = np.array([
        "|".join(z) for z in zip(out["read"], out["annotation"], out["sharing"])])
    return out

def adequacy(y, component):
    n, p, k = len(y), int(np.sum(y)), len(np.unique(component))
    return {"groups": n, "positive_groups": p, "components": k,
            "adequate": bool(n >= 200 and p >= 20 and k >= 20)}

def point_rows(d, strata):
    """All fixed model and stratum rows, including empty/single-class/inadequate cells."""
    y, comp, fold = (d[x].to_numpy() for x in ("target", "dependency_component_id", "inner_fold"))
    overall, foldrows, stratified = [], [], []
    scopes = {"ELIGIBLE_DEVELOPMENT": np.ones(len(d), dtype=bool), "EN71_COMMON_SUPPORT": d.EN71.notna().to_numpy()}
    for cohort, mask in scopes.items():
        for m in MODELS:
            valid = mask & np.isfinite(d[m].to_numpy())
            yv, sv, fv = y[valid], d[m].to_numpy()[valid], fold[valid]
            p, ff = Prepared(yv, sv, fv, m in PROBABILITIES).calculate(folds=True)
            head = {"cohort": cohort, "model": m, "eligible_groups": int(mask.sum()),
                    "scored_groups": int(valid.sum()),
                    "scoring_coverage": float(valid.sum()/mask.sum()) if mask.sum() else np.nan,
                    **adequacy(yv, comp[valid])}
            overall.append(dict(head, **p))
            for f, fm in enumerate(ff):
                scopefold = mask & (fold == f)
                validfold = valid & (fold == f)
                foldrows.append({"cohort": cohort, "model": m, "fold": f,
                    "eligible_groups": int(scopefold.sum()), **adequacy(y[validfold], comp[validfold]),
                    "scoring_coverage": float(validfold.sum()/scopefold.sum()) if scopefold.sum() else np.nan,
                    **fm})
    axes = dict(LEVELS, joint_read_annotation_sharing=JOINT)
    for axis, levels in axes.items():
        for level in levels:
            mask = strata[axis] == level
            for m in MODELS:
                valid = mask & np.isfinite(d[m].to_numpy())
                r = Prepared(y[valid], d[m].to_numpy()[valid], fold[valid], m in PROBABILITIES).calculate()
                stratified.append({"axis": axis, "level": level, "model": m,
                    "eligible_groups": int(mask.sum()), **adequacy(y[valid], comp[valid]),
                    "scoring_coverage": float(valid.sum()/mask.sum()) if mask.sum() else np.nan, **r})
    return overall, foldrows, stratified

def compare_prepare(d, pair, strata):
    a, b = pair
    ok = np.isfinite(d[a].to_numpy()) & np.isfinite(d[b].to_numpy())
    z = d.loc[ok]
    y, fold, comp = (z[c].to_numpy() for c in ("target", "inner_fold", "dependency_component_id"))
    _, ci = np.unique(comp, return_inverse=True)
    models = [Prepared(y, z[m].to_numpy(), fold, m in PROBABILITIES) for m in pair]
    joint = strata["joint_read_annotation_sharing"][ok]
    controls, rows = [], []
    for lev in JOINT:
        ix = np.flatnonzero(joint == lev)
        ad = adequacy(y[ix], comp[ix])
        aa, bb = [Ranked(y[ix], z[m].to_numpy()[ix], m in PROBABILITIES) for m in pair]
        ap, bp = aa.calculate(), bb.calculate()
        rows.append({"contrast": a + "_minus_" + b, "level": lev, **ad,
                     "a_enrichment": ap["full_ties_10_enrichment"],
                     "b_enrichment": bp["full_ties_10_enrichment"],
                     "delta": ap["full_ties_10_enrichment"] - bp["full_ties_10_enrichment"]})
        if ad["adequate"]:
            controls.append((len(ix), ix, aa, bb))
    return ok, z, ci, models, controls, rows

def controlled_delta(controls, w):
    den = sum(n for n, _, _, _ in controls)
    if not den:
        return np.nan
    values = [n * (a.calculate(w[ix])["full_ties_10_enrichment"] -
                   b.calculate(w[ix])["full_ties_10_enrichment"]) for n, ix, a, b in controls]
    # Never silently drop an undefined bootstrap cell or reweight to a different estimand.
    return float(np.sum(values) / den)

def contrast_metrics(a, b, controls, w):
    aa, bb = a.calculate(w), b.calculate(w)
    keys = [f"{mode}_{int(q*100):02d}_enrichment" for mode in MODES for q in Q]
    keys += ["AUC", "within_fold_pair_weighted_AUC"]
    if a.rank.probability and b.rank.probability:
        keys += ["Brier", "log_loss"]
    out = {k: aa[k] - bb[k] for k in keys}
    out["joint_n_weighted_full_ties_10_enrichment"] = controlled_delta(controls, w)
    return out

def interval(values):
    v = np.asarray(values, dtype=float)
    valid = v[np.isfinite(v)]
    ans = {"bootstrap_draws_requested": len(v), "bootstrap_draws_defined": len(valid),
           "bootstrap_draws_undefined": int(len(v)-len(valid)),
           "CI_status": "COMPLETE" if len(valid) == len(v) else "UNDEFINED_DRAWS_PRESENT"}
    # Incomplete draws are explicitly conditional-on-defined and may not support conclusions.
    for confidence, lo, hi in (("95", .025, .975), ("99_375", .003125, .996875)):
        ends = np.quantile(valid, [lo, hi], method="linear") if len(valid) else [np.nan, np.nan]
        ans[f"CI_{confidence}_low"], ans[f"CI_{confidence}_high"] = map(float, ends)
    return ans

def save_table(folder, name, rows):
    frame = pd.DataFrame(rows)
    p = folder / (name + ".tsv")
    with p.open("x") as f:
        frame.to_csv(f, sep="\t", index=False, na_rep="NA", float_format="%.17g")
    q = folder / (name + ".parquet")
    require(not q.exists(), "Refusing Parquet overwrite")
    with q.open("xb") as f:
        pq.write_table(pa.Table.from_pandas(frame, preserve_index=False), f, compression="zstd")
    return [p, q]

def evaluate(output_root):
    start, cpu0 = time.monotonic(), time.process_time()
    resource.setrlimit(resource.RLIMIT_CPU, (CPU_LIMIT, CPU_LIMIT))
    binding, input_sha = verify_bound_inputs()
    out = Path(output_root)
    require(out in {HERE / "results", HERE / "results_A", HERE / "results_B"}, "Only approved independent evaluation result roots are allowed")
    require(out.resolve().is_relative_to(HERE.resolve()), "Output path escape")
    require(not out.exists(), "Results path already exists: no overwrite/retry into accepted generation")
    # All outcomes are deserialized only after exact hashes and explicit CLI authorization pass.
    d, meta = load_data()
    strata = make_strata(meta)
    out.mkdir()
    paths = []
    overall, folds, stratified = point_rows(d, strata)
    for name, rows in (("model_metrics", overall), ("fold_metrics", folds), ("stratified_metrics", stratified)):
        paths += save_table(out, name, rows)
    contrasts, jointrows, commonfolds, commonmodels, commonstrata, bootstrap = [], [], [], [], [], []
    for pair in CONTRASTS:
        ok, z, ci, models, controls, rows = compare_prepare(d, pair, strata)
        jointrows.extend(rows)
        name = pair[0] + "_minus_" + pair[1]
        y, comp = z.target.to_numpy(), z.dependency_component_id.to_numpy()
        head = {"contrast": name, "a": pair[0], "b": pair[1], **adequacy(y, comp),
                "eligible_groups": len(d), "common_scoring_coverage": len(z)/len(d),
                "a_missing_groups": int(d[pair[0]].isna().sum()),
                "b_missing_groups": int(d[pair[1]].isna().sum()),
                "joint_fixed_adequate_cells": len(controls),
                "joint_fixed_weight_groups": sum(v[0] for v in controls)}
        w = np.ones(len(z))
        point = contrast_metrics(*models, controls, w)
        a_point, b_point = [m.calculate(w) for m in models]
        head.update(common_prevalence=float(np.mean(y)) if len(y) else np.nan,
                    common_group_ids_sha256=hashlib.sha256("\\n".join(sorted(z.index)).encode()).hexdigest(),
                    a_PPV=a_point["full_ties_10_PPV"], b_PPV=b_point["full_ties_10_PPV"],
                    a_actual_selection_coverage=a_point["full_ties_10_selection_coverage"],
                    b_actual_selection_coverage=b_point["full_ties_10_selection_coverage"])
        for axis, levels in LEVELS.items():
            for level in levels:
                ix = strata[axis][ok] == level
                for m in pair:
                    r = Prepared(y[ix], z[m].to_numpy()[ix], z.inner_fold.to_numpy()[ix],
                                 m in PROBABILITIES).calculate()
                    commonstrata.append({"contrast": name, "axis": axis, "level": level,
                                         "model": m, **adequacy(y[ix], comp[ix]), **r})
        for m, model in zip(pair, models):
            mp, ff = model.calculate(w, folds=True)
            commonmodels.append({"contrast": name, "model": m, **adequacy(y, comp), **mp})
            for f, fm in enumerate(ff):
                ix = z.inner_fold.eq(f).to_numpy()
                commonfolds.append({"contrast": name, "model": m, "fold": f,
                                    **adequacy(y[ix], comp[ix]), **fm})
        # Draw dependence components in THIS pair's common support; each component's
        # complete group multiplicity is retained. Same paired weights for both models.
        K = len(np.unique(ci))
        rng = np.random.default_rng(SEED)
        samples = {metric: [] for metric in point}
        for iteration in range(DRAWS):
            counts = np.bincount(rng.integers(K, size=K), minlength=K) if K else np.array([])
            weights = counts[ci].astype(float)
            delta = contrast_metrics(*models, controls, weights)
            for metric, value in delta.items():
                samples[metric].append(value)
            bootstrap.append({"contrast": name, "iteration": iteration, **delta})
            if iteration % 500 == 0:
                print(json.dumps({"event": "bootstrap", "contrast": name, "iteration": iteration,
                                  "cpu_seconds": time.process_time()-cpu0}), flush=True)
        for metric, value in point.items():
            contrasts.append({**head, "metric": metric, "delta_a_minus_b": value,
                "higher_delta_favors_a": metric not in ("Brier", "log_loss"),
                "family_size": 8, "family_policy": "Bonferroni within each metric; not across all sensitivities",
                "primary_metric": metric == "full_ties_10_enrichment",
                **interval(samples[metric])})
    for name, rows in (
        ("contrast_CI", contrasts), ("paired_common_support_models", commonmodels),
        ("paired_common_support_folds", commonfolds), ("paired_stratified_common_support", commonstrata),
        ("paired_joint_all_45_cells", jointrows),
        ("paired_cluster_bootstrap", bootstrap)):
        paths += save_table(out, name, rows)
    receipt = {
        "status": "COMPLETED_DESCRIPTIVE_DEVELOPMENT_EVALUATION",
        "evidence_status": EVIDENCE, "contract_sha256": CONTRACT_SHA,
        "input_binding_sha256": input_sha, "code_freeze_sha256": binding["code_freeze_sha256"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.monotonic()-start, "cpu_seconds": time.process_time()-cpu0,
        "bootstrap_draws": DRAWS, "bootstrap_seed": SEED,
        "intervals": "95% descriptive; 99.375% eight-contrast family adjusted; no confirmatory/equivalence inference",
        "independent_unit": "inherited dependence component; NOT independent biological replicate",
        "evaluation_unit": "one promoter group, one Bernoulli target",
        "log_loss_clip": LOG_EPS,
        "tie_semantics": "full exact ties primary; uniform fractional boundary at exact group budget secondary",
        "joint_weight_policy": "point-common-support n, fixed adequate cells n>=200/positive>=20/components>=20; all 45 cells retained",
        "new_fits": 0, "final_external_outcome_reads": 0,
        "undefined_estimates_serialized_as": "JSON null / TSV NA / Parquet IEEE NaN or null",
        "files": [file_record(str(p.relative_to(ROOT)), "DESCRIPTIVE_RESULT") for p in paths],
    }
    dump_new(out / "evaluation_complete.json", receipt)
    print(json.dumps({"status": receipt["status"], "cpu_seconds": receipt["cpu_seconds"]}))

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("interface", "bind-inputs", "evaluate"))
    p.add_argument("--output-root", default=str(HERE / "results"),
                   help="Only evaluation/results, evaluation/results_A or evaluation/results_B; must not exist")
    p.add_argument("--authorize-development-evaluation", action="store_true",
                   help="Explicit post-review authorization; required to deserialize development outcomes")
    args = p.parse_args()
    if args.command == "interface":
        print(json.dumps({"whitelist": whitelist(), "fixed_sha256": FIXED_HASHES,
                          "new_results_must_be_bound_before_evaluation": True,
                          "output": str(HERE / "results"), "contrasts": CONTRASTS}, indent=2))
    elif args.command == "bind-inputs":
        bind_inputs()
    else:
        require(args.authorize_development_evaluation, "No real evaluation without explicit reviewed authorization")
        evaluate(args.output_root)

if __name__ == "__main__":
    main()

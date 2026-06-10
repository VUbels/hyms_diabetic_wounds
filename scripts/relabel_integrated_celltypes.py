#!/usr/bin/env python3
"""
relabel_integrated_celltypes.py
###############################################################################
Re-label obs['predicted_cell_type'] on nichecompass_integrated.h5ad from the
per-sample  ./annotated_data/<sample>/<sample>_annotated.h5ad  files, WITHOUT
retraining NicheCompass.

Why this is safe (and why no retrain / no model touch is needed):
NicheCompass niches are a function of expression + the spatial graph + the gene-
program priors + the batch covariate. They never see the cell-type labels. So
adopting a new annotation is a pure column swap on the trained object; the niche
assignment, the latent embedding and the active gene programs are unchanged.
(See common_py_functions.relabel_from_source — this script is a thin, auditable
wrapper around it.)

Key facts wired in as defaults (override on the CLI if your setup differs):
  * The integrated cell-type column read by every niche_*/comm_* script is
    CELLTYPE_KEY = 'predicted_cell_type'  -> that is what we OVERWRITE.
  * The new (curated) labels live in the SOURCE files under 'cell_type'
    -> that is --label-col.
  * Integrated barcodes are assumed to be  '{sample_id}{sep}{barcode}', where
    sample_id is the per-sample folder name (e.g. XE765-L) and sep defaults to
    '_'. The 0.99 overlap guard refuses to write if that assumption is wrong,
    so a key-format mismatch fails loudly instead of corrupting labels.

Run from the PROJECT ROOT so the default relative paths resolve:
    python scripts/relabel_integrated_celltypes.py --check   # inspect, write nothing
    python scripts/relabel_integrated_celltypes.py           # apply (in place, with backup)
###############################################################################
"""

import os
import sys
import glob
import shutil
import argparse
import datetime

# Make common_py_functions importable no matter what the CWD is.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common_py_functions as cf          # noqa: E402
import pandas as pd                       # noqa: E402


def build_mapping(source, label_col, sep, barcode_col=None):
    """Reconstruct the {integrated_barcode: label} mapping for an annotated_dir,
    so we can report the overlap and the label shift BEFORE committing anything.

    The per-cell id is obs[barcode_col] when barcode_col is given, else the obs
    index. Keys are '{sample_id}{sep}{id}'."""
    files = (sorted(glob.glob(os.path.join(source, "*", "*_annotated.h5ad")))
             or sorted(glob.glob(os.path.join(source, "*_annotated.h5ad"))))
    if not files:
        raise FileNotFoundError(f"No *_annotated.h5ad found under {source}")
    mapping, per_file = {}, []
    for f in files:
        sid = os.path.basename(os.path.dirname(f))
        obs = cf._read_obs_only(f)                      # cheap, backed read of .obs only
        if label_col not in obs.columns:
            per_file.append((sid, f, 0,
                             f"MISSING label '{label_col}'. obs columns: {list(obs.columns)}"))
            continue
        if barcode_col is not None and barcode_col not in obs.columns:
            per_file.append((sid, f, 0,
                             f"MISSING barcode '{barcode_col}'. obs columns: {list(obs.columns)}"))
            continue
        ids = (obs[barcode_col] if barcode_col is not None else obs.index).astype(str)
        n = 0
        for bc, lab in zip(ids, obs[label_col].astype(str)):
            mapping[f"{sid}{sep}{bc}"] = lab
            n += 1
        per_file.append((sid, f, n, "ok"))
    return mapping, files, per_file


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--h5ad", default=cf.DEFAULT_H5AD,
                    help=f"Integrated object to update (default: {cf.DEFAULT_H5AD}).")
    ap.add_argument("--source", default="./annotated_data",
                    help="Dir of per-sample <sid>/<sid>_annotated.h5ad (default: ./annotated_data).")
    ap.add_argument("--label-col", "--label_col", dest="label_col", default="cell_type",
                    help="Cell-type column IN THE SOURCE files (default: cell_type).")
    ap.add_argument("--celltype-key", "--celltype_key", dest="celltype_key",
                    default=cf.CELLTYPE_KEY,
                    help=f"Column to OVERWRITE in the integrated object (default: {cf.CELLTYPE_KEY}).")
    ap.add_argument("--sep", default="_",
                    help="Separator in integrated obs_names '{sample_id}{sep}{barcode}' (default: '_').")
    ap.add_argument("--barcode-col", "--barcode_col", dest="barcode_col", default=None,
                    help="Per-sample obs column to use as the barcode instead of the obs index "
                         "(set this if the JOIN-KEY test shows e.g. 'original_cell_id' matches, "
                         "not the index).")
    ap.add_argument("--min-overlap", "--min_overlap", dest="min_overlap",
                    type=float, default=0.99,
                    help="Refuse to write below this matched fraction (default: 0.99).")
    ap.add_argument("--out", default=None,
                    help="Write here instead of in place. Default: overwrite --h5ad.")
    ap.add_argument("--check", action="store_true",
                    help="Dry run: report overlap + label shift, write nothing.")
    ap.add_argument("--no-backup", action="store_true",
                    help="Skip the timestamped .bak copy before an in-place write.")
    args = ap.parse_args()

    import anndata as ad                                # lazy so --help works without scanpy stack

    h5ad = cf.resolve_input_path(args.h5ad, what="integrated .h5ad")
    source = cf.resolve_input_path(args.source, what="annotated_data dir")
    out = args.out or h5ad
    in_place = os.path.abspath(out) == os.path.abspath(h5ad)

    print(f"[RELABEL] integrated : {h5ad}")
    print(f"[RELABEL] source dir : {source}")
    print(f"[RELABEL] mapping    : source['{args.label_col}']  ->  integrated['{args.celltype_key}']")
    print(f"[RELABEL] key format : '{{sample_id}}{args.sep}{{{'barcode' if args.barcode_col is None else args.barcode_col}}}'  (sample_id = folder name)")

    # ---------- pre-flight (no writes) ----------
    mapping, files, per_file = build_mapping(source, args.label_col, args.sep, args.barcode_col)
    print(f"\n[PREFLIGHT] {len(files)} annotated file(s):")
    for sid, f, n, status in per_file:
        print(f"    {sid:<18}{n:>9} cells  [{status}]  {os.path.relpath(f)}")

    adata = ad.read_h5ad(h5ad)
    names = pd.Series(adata.obs_names.astype(str), index=adata.obs_names)   # Series so .notna().values is safe
    new = names.map(mapping)
    matched = new.notna()
    overlap = float(matched.mean())

    print(f"\n[PREFLIGHT] integrated cells        : {adata.n_obs}")
    print(f"[PREFLIGHT] matched by key           : {int(matched.sum())} ({overlap:.2%})")
    print(f"[PREFLIGHT] proposed label classes   : {int(new.dropna().nunique())}")
    print(f"[PREFLIGHT] example integrated names : {list(names.values[:3])}")
    print(f"[PREFLIGHT] example source keys      : {list(mapping)[:3]}")
    if (~matched).any():
        print(f"[PREFLIGHT] example UNMATCHED names  : {list(names.values[(~matched).values][:5])}")

    if overlap < args.min_overlap:
        print(f"\n[ABORT] overlap {overlap:.2%} < --min-overlap {args.min_overlap:.0%}. Nothing written.\n"
              f"        Compare the example integrated names vs example source keys above:\n"
              f"        adjust --label-col / --sep / --barcode-col so the keys line up\n"
              f"        (the diagnostic JOIN-KEY test shows which source id matches).")
        sys.exit(2)

    # label-shift crosstab (only reached when overlap is sufficient -> matched > 0)
    if args.celltype_key in adata.obs.columns:
        old = adata.obs[args.celltype_key].astype(str)
        m = matched.values
        moved = int((old.values[m] != new.astype(str).values[m]).sum())
        denom = max(int(m.sum()), 1)
        print(f"\n[PREFLIGHT] labels that CHANGE       : {moved} of {int(m.sum())} matched ({moved/denom:.1%})")
        ct = pd.crosstab(old[matched], new[matched].astype(str))
        with pd.option_context("display.max_rows", 60, "display.max_columns", 60,
                               "display.width", 220):
            print("[PREFLIGHT] old (rows) x new (cols) crosstab:")
            print(ct)
    else:
        print(f"\n[PREFLIGHT] integrated has no '{args.celltype_key}' yet; it will be created.")

    if args.check:
        print("\n[CHECK] dry run only — nothing written. Re-run without --check to apply.")
        return

    # ---------- apply (reuse the validated mapping; supports --barcode-col) ----------
    if args.celltype_key in adata.obs.columns:
        adata.obs[args.celltype_key + "_prev"] = adata.obs[args.celltype_key].astype(str).values
    adata.obs[args.celltype_key] = pd.Categorical(new.astype(str).values)

    tmp = out + ".tmp_relabel"
    adata.write_h5ad(tmp)
    if in_place and not args.no_backup:
        stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        bak = f"{h5ad}.bak-{stamp}"
        shutil.copy2(h5ad, bak)
        print(f"[BACKUP] {bak}")
    os.replace(tmp, out)
    print(f"[WROTE]  {out}")
    print(f"[DONE]   obs['{args.celltype_key}'] updated ({int(new.dropna().nunique())} classes); "
          f"previous labels preserved in obs['{args.celltype_key}_prev'].")


if __name__ == "__main__":
    main()

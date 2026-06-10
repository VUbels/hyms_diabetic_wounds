#!/usr/bin/env python
"""
set_condition.py
###############################################################################
STEP 0 of the pipeline. Run this once, before any analysis script.

It writes a single standardised condition column on the integrated object whose
only two values are the literal strings "CNTRL" and "TEST". Every downstream
script then reads that column and reports every contrast as TEST vs CNTRL, so
there is never a per-script --cntrl/--test flag to keep consistent.

You define the mapping in one of two ways:

  (A) from an existing obs column
      python set_condition.py \
          --h5ad nichecompass_results/objects/nichecompass_integrated.h5ad \
          --from-key sample_id \
          --cntrl-value XE765-L XE766-L XE768-L \
          --test-value XE789-D XE791-D XE793-D

  (B) by listing the sample IDs in each arm
      python set_condition.py \
          --h5ad nichecompass_results/objects/nichecompass_integrated.h5ad \
          --cntrl-samples XE765-L XE766-L XE768-L \
          --test-samples  XE789-D XE791-D XE793-D

    Categories (6, object): ['XE765-L', 'XE766-L', 'XE768-L', 'XE789-D', 'XE791-D', 'XE793-D']
    
By default the object is written back in place (the path the analysis scripts
read), and any prior values in the condition column are preserved in a
'<condition>_original' column. Cells that map to neither arm are left as NaN and
reported; pass --drop-unmapped to remove them instead.
###############################################################################
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd

import common_py_functions as cf


###############################################################################
# ARGUMENTS
###############################################################################
def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h5ad", default=cf.DEFAULT_H5AD, metavar="FILE",
                   help=f"Integrated .h5ad to annotate (default: {cf.DEFAULT_H5AD}). "
                        "This is the consistent output of the integration step.")
    p.add_argument("--out", default=None, metavar="FILE",
                   help="Where to write the annotated object (default: in place, "
                        "overwriting --h5ad). Downstream scripts read this object.")
    p.add_argument("--condition-key", "--condition_key", dest="condition_key",
                   default=cf.CONDITION_KEY,
                   help=f"obs column to write (default: {cf.CONDITION_KEY}).")
    p.add_argument("--sample-key", "--sample_key", dest="sample_key",
                   default=cf.SAMPLE_KEY,
                   help=f"obs column with sample IDs (default: {cf.SAMPLE_KEY}).")

    # mode A: map from an existing obs column
    p.add_argument("--from-key", "--from_key", dest="from_key", default=None,
                   help="Existing obs column to map from (mode A).")
    p.add_argument("--cntrl-value", "--cntrl_value", dest="cntrl_value",
                   nargs="*", default=None,
                   help="Value(s) of --from-key that become CNTRL (mode A).")
    p.add_argument("--test-value", "--test_value", dest="test_value",
                   nargs="*", default=None,
                   help="Value(s) of --from-key that become TEST (mode A).")

    # mode B: assign by sample membership
    p.add_argument("--cntrl-samples", "--cntrl_samples", dest="cntrl_samples",
                   nargs="*", default=None,
                   help="Sample IDs that become CNTRL (mode B).")
    p.add_argument("--test-samples", "--test_samples", dest="test_samples",
                   nargs="*", default=None,
                   help="Sample IDs that become TEST (mode B).")

    p.add_argument("--drop-unmapped", "--drop_unmapped", dest="drop_unmapped",
                   action="store_true",
                   help="Drop cells that map to neither arm (default: keep as NaN).")
    return p.parse_args()


###############################################################################
# MAPPING
###############################################################################
def build_condition(adata, args):
    """Return a pandas Series (index = obs_names) of CNTRL / TEST / NaN."""
    cntrl, test = cf.CONDITION_CNTRL, cf.CONDITION_TEST
    n = adata.n_obs
    out = pd.Series([np.nan] * n, index=adata.obs_names, dtype=object)

    mode_a = args.from_key is not None
    mode_b = args.cntrl_samples is not None or args.test_samples is not None
    if mode_a and mode_b:
        raise SystemExit("Choose ONE mapping mode: --from-key (A) OR "
                         "--cntrl-samples/--test-samples (B), not both.")
    if not mode_a and not mode_b:
        raise SystemExit("No mapping given. Use --from-key + --cntrl-value/"
                         "--test-value (A), or --cntrl-samples/--test-samples (B).")

    if mode_a:
        if not args.cntrl_value or not args.test_value:
            raise SystemExit("Mode A needs both --cntrl-value and --test-value.")
        if args.from_key not in adata.obs.columns:
            raise SystemExit(
                f"--from-key '{args.from_key}' not in obs. Available: "
                f"{list(adata.obs.columns)[:20]}")
        col = adata.obs[args.from_key].astype(str)
        cset = set(map(str, args.cntrl_value))
        tset = set(map(str, args.test_value))
        overlap = cset & tset
        if overlap:
            raise SystemExit(f"Values appear in both arms: {sorted(overlap)}")
        out[col.isin(cset).values] = cntrl
        out[col.isin(tset).values] = test
        present = set(col.unique())
        missing = (cset | tset) - present
        if missing:
            print(f"  [WARN] requested values not found in '{args.from_key}': "
                  f"{sorted(missing)}")
    else:
        if args.sample_key not in adata.obs.columns:
            raise SystemExit(
                f"--sample-key '{args.sample_key}' not in obs. Available: "
                f"{list(adata.obs.columns)[:20]}")
        samp = adata.obs[args.sample_key].astype(str)
        cset = set(map(str, args.cntrl_samples or []))
        tset = set(map(str, args.test_samples or []))
        overlap = cset & tset
        if overlap:
            raise SystemExit(f"Samples listed in both arms: {sorted(overlap)}")
        out[samp.isin(cset).values] = cntrl
        out[samp.isin(tset).values] = test
        present = set(samp.unique())
        missing = (cset | tset) - present
        if missing:
            print(f"  [WARN] requested sample IDs not present: {sorted(missing)}")
    return out


###############################################################################
# MAIN
###############################################################################
def main():
    args = parse_args()
    import anndata as ad

    path = cf.resolve_input_path(args.h5ad, required=True, what="--h5ad")
    out_path = args.out or path
    cf.banner("SET CONDITION (CNTRL vs TEST)")
    print(f"[LOAD] {path}")
    adata = ad.read_h5ad(path)
    print(f"  {adata.n_obs} cells")

    cond = build_condition(adata, args)

    n_cntrl = int((cond == cf.CONDITION_CNTRL).sum())
    n_test = int((cond == cf.CONDITION_TEST).sum())
    n_unmapped = int(cond.isna().sum())
    print(f"  CNTRL cells: {n_cntrl}")
    print(f"  TEST  cells: {n_test}")
    print(f"  unmapped   : {n_unmapped}")
    if n_cntrl == 0 or n_test == 0:
        raise SystemExit("One arm is empty after mapping - check your values/IDs.")

    # preserve any prior values in the target column
    key = args.condition_key
    if key in adata.obs.columns:
        backup = f"{key}_original"
        if backup not in adata.obs.columns:
            adata.obs[backup] = adata.obs[key].astype(str).values
            print(f"  prior '{key}' preserved as '{backup}'")

    adata.obs[key] = cond.values

    if n_unmapped:
        if args.drop_unmapped:
            keep = adata.obs[key].notna().values
            adata = adata[keep].copy()
            print(f"  dropped {n_unmapped} unmapped cells -> {adata.n_obs} remain")
        else:
            print(f"  [WARN] {n_unmapped} cells left as NaN in '{key}' "
                  f"(use --drop-unmapped to remove them)")

    # sample-level summary (the unit of analysis)
    smap = (adata.obs.groupby(args.sample_key)[key]
            .agg(lambda s: s.dropna().iloc[0] if s.notna().any() else "NaN"))
    print("  sample -> condition:")
    for s, c in smap.items():
        print(f"    {s}: {c}")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    adata.write_h5ad(out_path)
    print(f"[WRITE] {out_path}")
    print(f"  obs['{key}'] now holds '{cf.CONDITION_CNTRL}' / '{cf.CONDITION_TEST}'. "
          f"Downstream scripts need no --cntrl/--test.")


if __name__ == "__main__":
    main()

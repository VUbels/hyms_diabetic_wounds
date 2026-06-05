#!/usr/bin/env python
"""
Rename condition labels in the NicheCompass-integrated AnnData object.

    healthy  -> high_senescence
    chronic  -> low_senescence

The file is overwritten in place at:
    nichecompass_results/objects/nichecompass_integrated.h5ad
"""

from pathlib import Path
import anndata as ad

H5AD_PATH = Path("nichecompass_results/objects/nichecompass_integrated.h5ad")
CONDITION_COL = "condition"  # <-- change if your obs column has a different name
RENAME_MAP = {"healthy": "high_senescence", "chronic": "low_senescence"}


def main() -> None:
    print(f"[LOAD] {H5AD_PATH}")
    adata = ad.read_h5ad(H5AD_PATH)

    if CONDITION_COL not in adata.obs.columns:
        raise KeyError(
            f"Column '{CONDITION_COL}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    print(f"  Before: {adata.obs[CONDITION_COL].value_counts().to_dict()}")

    # Works whether the column is categorical or plain object/string.
    series = adata.obs[CONDITION_COL]
    if hasattr(series, "cat"):
        # Preserve categorical dtype but rename categories cleanly.
        # Using map + astype('category') handles cases where new names
        # might collide or where categories need to be dropped.
        new_values = series.astype(str).replace(RENAME_MAP)
        adata.obs[CONDITION_COL] = new_values.astype("category")
    else:
        adata.obs[CONDITION_COL] = series.replace(RENAME_MAP)

    print(f"  After:  {adata.obs[CONDITION_COL].value_counts().to_dict()}")

    print(f"[SAVE] {H5AD_PATH}")
    adata.write_h5ad(H5AD_PATH, compression="gzip")
    print("[DONE]")


if __name__ == "__main__":
    main()
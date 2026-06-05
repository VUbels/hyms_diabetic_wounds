import gzip, json
import pandas as pd
import numpy as np
from pathlib import Path
from shapely.geometry import shape, Polygon
from scipy.spatial import cKDTree

SAMPLE = "XE765-L"
XEN_DIR = Path(f"/mnt/d/hyms/metal_diabetes/{SAMPLE}")
PRO_DIR = Path(f"./proseg_results_path3/{SAMPLE}")

# ── 1. Inspect proseg's cell-metadata ─────────────────────────────────
pro_meta = pd.read_csv(PRO_DIR / "cell-metadata.csv.gz")
print(f"Proseg cell-metadata shape: {pro_meta.shape}")
print(f"Proseg cell-metadata columns: {list(pro_meta.columns)}")
print(pro_meta.head(3))
print()

# ── 2. Xenium centroids and reconstructed polygons ────────────────────
xen_cells = pd.read_parquet(XEN_DIR / "cells.parquet")
xen_cells["cell_id"] = xen_cells["cell_id"].astype(str)
xen_centroids = xen_cells.set_index("cell_id")[["x_centroid", "y_centroid"]]

bounds = pd.read_parquet(XEN_DIR / "cell_boundaries.parquet")
bounds["cell_id"] = bounds["cell_id"].astype(str)
xenium_geoms = {}
for cid, grp in bounds.groupby("cell_id", sort=False):
    pts = list(zip(grp["vertex_x"].values, grp["vertex_y"].values))
    if len(pts) >= 3:
        p = Polygon(pts)
        if not p.is_valid:
            p = p.buffer(0)
        if p.is_valid and not p.is_empty:
            xenium_geoms[cid] = p
print(f"Xenium polygons: {len(xenium_geoms):,}")

# ── 3. Proseg polygons (keyed by integer `cell`) ──────────────────────
with gzip.open(PRO_DIR / "cell-polygons.geojson.gz", "rt") as f:
    gj = json.load(f)
proseg_geoms_by_int = {}
for feat in gj["features"]:
    cid = int(feat["properties"]["cell"])
    g = shape(feat["geometry"])
    if not g.is_valid:
        g = g.buffer(0)
    if g.is_valid and not g.is_empty:
        proseg_geoms_by_int[cid] = g
print(f"Proseg polygons: {len(proseg_geoms_by_int):,}")
print(f"Proseg cell integer range: {min(proseg_geoms_by_int)}..{max(proseg_geoms_by_int)}")
print()

# ── 4. Build proseg_int → xenium_cell_id mapping ──────────────────────
# Strategy A: direct mapping if cell-metadata has both columns
mapping = None
str_cols = [c for c in pro_meta.columns
            if pro_meta[c].dtype == object and c not in ("color",)]
for c in str_cols:
    sample_vals = pro_meta[c].dropna().head(5).astype(str).tolist()
    looks_like_xenium = any("-" in v and len(v) > 5 for v in sample_vals)
    if looks_like_xenium:
        print(f"Found xenium-style IDs in column '{c}': {sample_vals[:3]}")
        int_col = "cell" if "cell" in pro_meta.columns else pro_meta.columns[0]
        mapping = dict(zip(pro_meta[int_col].astype(int),
                           pro_meta[c].astype(str)))
        print(f"Direct mapping built from cell-metadata: {len(mapping):,} entries")
        break

# Strategy B: centroid matching (fallback)
if mapping is None:
    print("No direct ID column found; using centroid matching.")
    # Proseg centroid column names vary; detect them
    cand_x = [c for c in pro_meta.columns if "centroid_x" in c.lower() or c.lower() in ("x", "centroid_x")]
    cand_y = [c for c in pro_meta.columns if "centroid_y" in c.lower() or c.lower() in ("y", "centroid_y")]
    print(f"Proseg centroid columns detected: x={cand_x}, y={cand_y}")
    if not cand_x or not cand_y:
        raise SystemExit("Cannot find proseg centroid columns — inspect metadata above.")
    pro_xy = pro_meta[[cand_x[0], cand_y[0]]].values
    int_col = "cell" if "cell" in pro_meta.columns else pro_meta.columns[0]
    pro_ints = pro_meta[int_col].astype(int).values

    xen_xy = xen_centroids[["x_centroid", "y_centroid"]].values
    xen_ids = xen_centroids.index.values

    tree = cKDTree(xen_xy)
    dists, idx = tree.query(pro_xy, k=1)
    print(f"Centroid-match distance distribution: "
          f"median {np.median(dists):.2f} µm, 95th {np.percentile(dists, 95):.2f} µm, "
          f"max {dists.max():.2f} µm")
    valid = dists < 5.0  # 5 µm tolerance
    print(f"Matches within 5 µm: {valid.sum():,} / {len(pro_ints):,}")
    mapping = dict(zip(pro_ints[valid], xen_ids[idx[valid]]))
    print(f"Centroid-based mapping built: {len(mapping):,} entries")

# ── 5. Compute IoU on mapped pairs ────────────────────────────────────
rows = []
for pint, xcid in mapping.items():
    if pint not in proseg_geoms_by_int or xcid not in xenium_geoms:
        continue
    xg, pg = xenium_geoms[xcid], proseg_geoms_by_int[pint]
    u = xg.union(pg).area
    if u == 0:
        continue
    rows.append({
        "iou": xg.intersection(pg).area / u,
        "area_xen": xg.area,
        "area_pro": pg.area,
        "area_ratio": pg.area / xg.area if xg.area > 0 else np.nan,
    })
res = pd.DataFrame(rows)
print(f"\nIoU pairs computed: {len(res):,}")
print("\nIoU distribution (Xenium vs ProSeg, per matched cell):")
print(res["iou"].describe(percentiles=[0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95]))
print("\nArea ratio (ProSeg / Xenium):")
print(res["area_ratio"].describe(percentiles=[0.10, 0.50, 0.90]))
print(f"\nMean IoU: {res['iou'].mean():.3f}")
print(f"Fraction IoU > 0.5: {(res['iou'] > 0.5).mean():.1%}")
print(f"Fraction IoU > 0.8: {(res['iou'] > 0.8).mean():.1%}")

# Save the mapping for later use in hybrid construction
if mapping:
    pd.DataFrame({
        "proseg_int_id":   list(mapping.keys()),
        "xenium_cell_id":  list(mapping.values()),
    }).to_csv(PRO_DIR / "proseg_to_xenium_id_map.csv", index=False)
    print(f"\nMapping saved to {PRO_DIR / 'proseg_to_xenium_id_map.csv'}")
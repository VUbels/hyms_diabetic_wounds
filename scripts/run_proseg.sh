#!/usr/bin/env bash
#
# run_proseg.sh — Batch proseg pipeline for Xenium spatial transcriptomics
#
# Recursively finds all transcripts.parquet files under a target directory,
# runs proseg on each, and produces outputs compatible with:
#   - Python (spatialdata zarr)
#   - R / Seurat v5 (mtx + csv.gz metadata)
#   - Xenium Explorer (proseg-to-baysor conversion)
#
# Usage:
#   ./run_proseg.sh \
#     -i ./xenium_data \
#     -o ./proseg_results \
#     --ncomponents 15 \
#     --voxel-size 0.5 \
#     --cell-compactness 0.05 \
#     -t 32
#
# Options:
#   -i,  --input DIR                  Root directory containing Xenium region folders
#                                     (default: current working directory)
#   -o,  --output DIR                 Output root directory
#                                     (default: current working directory)
#   -t,  --threads N                  Number of threads per proseg run
#                                     (default: all available)
#   -p,  --parallel N                 Run N regions in parallel
#                                     (default: 1 = sequential)
#   -n,  --dry-run                    Print what would be run without executing
#        --proseg PATH                Path to proseg binary
#                                     (default: searches PATH, then
#                                      ./proseg/target/release/proseg)
#
# ── General / sampling arguments ────────────────────────────────────────────
#        --voxel-size S               Voxel size in microns on x/y axis
#                                     (proseg default: 0.5)
#        --burnin-voxel-size S        Larger voxel size for burn-in phase;
#                                     must be integer multiple of --voxel-size
#        --voxel-layers N             Number of z-axis voxel layers
#                                     (proseg default: 1)
#        --samples N                  Number of sampling iterations
#        --burnin-samples N           Number of burn-in sampling iterations
#
# ── Output arguments ────────────────────────────────────────────────────────
#        --output-expected-counts     Also write expected (non-integer) count
#                                     matrix to expected-counts.mtx.gz
#        --output-rates               Also write cell-by-gene Poisson rate
#                                     parameters to rates.csv.gz
#        --output-polygon-layers      Also write per-z-layer non-overlapping
#                                     polygons to cell-polygons-layers.geojson.gz
#        --exclude-transcripts        Exclude transcript positions from the
#                                     spatialdata zarr (saves disk space)
#        --compare                    Compare original vs proseg segmentation
#                                     after each region (default: on)
#        --no-compare                 Skip segmentation comparison
#
# ── Model arguments ─────────────────────────────────────────────────────────
#        --ncomponents N              NB mixture components for gene expression
#                                     (proseg default: 10)
#        --no-diffusion               Disable RNA diffusion/leakage modeling
#        --diffusion-probability P    Prior prob. transcript is diffused
#                                     (proseg default: 0.2)
#        --diffusion-sigma-far S      Prior SD on repositioning of diffused
#                                     transcripts (proseg default: 4)
#        --diffusion-sigma-near S     Prior SD on repositioning of
#                                     non-diffused transcripts (proseg default: 1)
#        --nuclear-reassignment-prob P  Prior prob. initial nuclear assignment
#                                       is incorrect (proseg default: 0.2)
#        --cell-compactness V         Voxel morphology compactness prior;
#                                     larger = less spherical cells
#                                     (proseg default: 0.03)
#
# ── Xenium-mask prior mode (alternative to transcript-table prior) ──────────
#        --use-xenium-mask            Rasterise Xenium cell_boundaries.parquet
#                                     into a dense uint32 mask and feed it to
#                                     proseg via --cellpose-masks. This makes
#                                     proseg's prior cover the full Xenium
#                                     polygon footprint (including anucleate
#                                     cells), unlike the default transcript-
#                                     table prior which only covers cells with
#                                     ≥1 assigned transcript.
#        --mask-scale S               Microns per pixel for the rasterised mask
#                                     (default: 0.5; typically matches
#                                     --voxel-size)
#
#   -h,  --help                       Show this help message
#
# Requirements:
#   - proseg v3+ (built via `cargo install proseg` or from source)
#   - proseg-to-baysor (installed alongside proseg)
#   - Python 3 + pandas + pyarrow (optional, for --compare)
#   - matplotlib (optional, for visual border comparison PDFs)
#
# Output structure per region:
#   <output_dir>/<region_name>/
#     ├── proseg-output.zarr                    # spatialdata format
#     ├── counts.mtx.gz                         # integer count matrix
#     ├── expected-counts.mtx.gz                # (optional) expected counts
#     ├── cell-metadata.csv.gz
#     ├── gene-metadata.csv.gz
#     ├── transcript-metadata.csv.gz
#     ├── rates.csv.gz                          # (optional)
#     ├── cell-polygons.geojson.gz              # 2D consensus polygons
#     ├── cell-polygons-layers.geojson.gz       # (optional) per-z-layer polygons
#     ├── segmentation_comparison.tsv            # (optional) original vs proseg stats
#     ├── segmentation_comparison_detail_200um.pdf    # (optional) zoomed border comparison
#     ├── segmentation_comparison_neighbourhood_500um.pdf  # (optional) wider border comparison
#     ├── segmentation_comparison_distributions.pdf   # (optional) area + transcript histograms
#     ├── proseg.log
#     └── xenium_explorer/
#         ├── transcript-metadata.csv
#         ├── cell-polygons.geojson
#         └── import_to_xenium_ranger.sh
#

set -euo pipefail

# ============================================================================
# Defaults
# ============================================================================
INPUT_DIR="$(pwd)"
OUTPUT_DIR="$(pwd)"
THREADS=""
PARALLEL_JOBS=1
DRY_RUN=false
PROSEG_BIN=""

# Flags that control which optional output files are requested
OUTPUT_EXPECTED_COUNTS=false
OUTPUT_RATES=false
OUTPUT_POLYGON_LAYERS=false
EXCLUDE_TRANSCRIPTS=false
COMPARE_SEG=true

# All remaining proseg model/sampling args are accumulated here and
# forwarded verbatim to the proseg binary.
EXTRA_PROSEG_ARGS=()

# Cellpose-mask mode: rasterize Xenium polygons and feed as dense prior
USE_XENIUM_MASK=false
MASK_SCALE=0.5

# ============================================================================
# Color output helpers
# ============================================================================
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()  { echo -e "${BLUE}[INFO]${NC}  $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC}    $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC}  $(date '+%Y-%m-%d %H:%M:%S') $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $(date '+%Y-%m-%d %H:%M:%S') $*"; }

# ============================================================================
# Argument parsing
# ============================================================================
show_help() {
    sed -n '/^# Usage:/,/^# Requirements:/{ /^# Requirements:/d; s/^# \?//p }' "$0"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in

        # ── Script-level arguments ───────────────────────────────────────────
        -i|--input)         INPUT_DIR="$2";   shift 2 ;;
        -o|--output)        OUTPUT_DIR="$2";  shift 2 ;;
        -t|--threads)       THREADS="$2";     shift 2 ;;
        -p|--parallel)      PARALLEL_JOBS="$2"; shift 2 ;;
        -n|--dry-run)       DRY_RUN=true;     shift   ;;
        --proseg)           PROSEG_BIN="$2";  shift 2 ;;
        -h|--help)          show_help ;;

        # ── Output flags (script intercepts these to build correct paths) ────
        --output-expected-counts)   OUTPUT_EXPECTED_COUNTS=true; shift ;;
        --output-rates)             OUTPUT_RATES=true;           shift ;;
        --output-polygon-layers)    OUTPUT_POLYGON_LAYERS=true;  shift ;;
        --exclude-transcripts)      EXCLUDE_TRANSCRIPTS=true;    shift ;;
        --compare)                  COMPARE_SEG=true;            shift ;;
        --no-compare)               COMPARE_SEG=false;           shift ;;

        # ── General / sampling arguments ────────────────
        --voxel-size)
            EXTRA_PROSEG_ARGS+=("--voxel-size" "$2"); shift 2 ;;
        --burnin-voxel-size)
            EXTRA_PROSEG_ARGS+=("--burnin-voxel-size" "$2"); shift 2 ;;
        --voxel-layers)
            EXTRA_PROSEG_ARGS+=("--voxel-layers" "$2"); shift 2 ;;
        --samples)
            EXTRA_PROSEG_ARGS+=("--samples" "$2"); shift 2 ;;
        --burnin-samples)
            EXTRA_PROSEG_ARGS+=("--burnin-samples" "$2"); shift 2 ;;

        # ── Model arguments ─────────────────────────────
        --ncomponents)
            EXTRA_PROSEG_ARGS+=("--ncomponents" "$2"); shift 2 ;;
        --prior-seg-reassignment-prob)
            EXTRA_PROSEG_ARGS+=("--prior-seg-reassignment-prob" "$2"); shift 2 ;;
        --max-transcript-nucleus-distance)
            EXTRA_PROSEG_ARGS+=("--max-transcript-nucleus-distance" "$2"); shift 2 ;;
        --enforce-connectivity)
            EXTRA_PROSEG_ARGS+=("--enforce-connectivity"); shift ;;
        --use-xenium-mask)
            USE_XENIUM_MASK=true; shift ;;
        --mask-scale)
            MASK_SCALE="$2"; shift 2 ;;
        --min-qv)
            EXTRA_PROSEG_ARGS+=("--min-qv" "$2"); shift 2 ;;
        --no-diffusion)
            EXTRA_PROSEG_ARGS+=("--no-diffusion"); shift ;;
        --diffusion-probability)
            EXTRA_PROSEG_ARGS+=("--diffusion-probability" "$2"); shift 2 ;;
        --diffusion-sigma-far)
            EXTRA_PROSEG_ARGS+=("--diffusion-sigma-far" "$2"); shift 2 ;;
        --diffusion-sigma-near)
            EXTRA_PROSEG_ARGS+=("--diffusion-sigma-near" "$2"); shift 2 ;;
        --nuclear-reassignment-prob)
            # proseg uses underscore in the flag name
            EXTRA_PROSEG_ARGS+=("--nuclear-reassignment-prob" "$2"); shift 2 ;;
        --cell-compactness)
            EXTRA_PROSEG_ARGS+=("--cell-compactness" "$2"); shift 2 ;;

        *)
            log_error "Unknown argument: $1"
            echo "       Run with -h / --help to see all options."
            exit 1 ;;
    esac
done

# ============================================================================
# Resolve proseg binary
# ============================================================================
resolve_proseg() {
    if [[ -n "$PROSEG_BIN" ]]; then
        if [[ -x "$PROSEG_BIN" ]]; then
            echo "$PROSEG_BIN"
            return
        else
            log_error "Specified proseg binary not found or not executable: $PROSEG_BIN"
            exit 1
        fi
    fi

    if command -v proseg &>/dev/null; then
        command -v proseg
        return
    fi

    local local_bin="./proseg/target/release/proseg"
    if [[ -x "$local_bin" ]]; then
        echo "$local_bin"
        return
    fi

    log_error "proseg binary not found. Install with 'cargo install proseg' or build from source."
    exit 1
}

PROSEG_BIN="$(resolve_proseg)"
log_info "Using proseg binary: $PROSEG_BIN"

PROSEG_TO_BAYSOR=""
if command -v proseg-to-baysor &>/dev/null; then
    PROSEG_TO_BAYSOR="$(command -v proseg-to-baysor)"
elif [[ -x "./proseg/target/release/proseg-to-baysor" ]]; then
    PROSEG_TO_BAYSOR="./proseg/target/release/proseg-to-baysor"
else
    log_warn "proseg-to-baysor not found. Xenium Explorer output will be skipped."
fi

# ============================================================================
# Validate input directory
# ============================================================================
if [[ ! -d "$INPUT_DIR" ]]; then
    log_error "Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# ============================================================================
# Discover all transcripts.parquet files
# ============================================================================
log_info "Scanning for transcripts.parquet files in: $INPUT_DIR"

mapfile -t PARQUET_FILES < <(find -L "$INPUT_DIR" -name "transcripts.parquet" -type f | sort)

if [[ ${#PARQUET_FILES[@]} -eq 0 ]]; then
    log_error "No transcripts.parquet files found under $INPUT_DIR"
    exit 1
fi

log_info "Found ${#PARQUET_FILES[@]} region(s):"
for f in "${PARQUET_FILES[@]}"; do
    region_dir="$(dirname "$f")"
    region_name="$(basename "$region_dir")"
    echo "         - $region_name"
done

# ============================================================================
# Compare original Xenium segmentation vs proseg resegmentation
#
# Reads the original Xenium cell/transcript data and proseg outputs,
# then prints a structured comparison, writes a TSV summary, and
# generates a multi-panel PDF comparing polygon borders visually.
#
# Requires Python 3 with pyarrow, pandas, numpy.
# Visual comparison additionally requires matplotlib.
# ============================================================================
compare_segmentation() {
    local region_dir="$1"    # original Xenium region directory
    local out_base="$2"      # proseg output directory for this region
    local region_name="$3"

    local xenium_transcripts="$region_dir/transcripts.parquet"
    local xenium_cells="$region_dir/cells.parquet"
    local xenium_boundaries="$region_dir/cell_boundaries.parquet"
    local proseg_cell_meta="$out_base/cell-metadata.csv.gz"
    local proseg_tx_meta="$out_base/transcript-metadata.csv.gz"
    local proseg_polygons="$out_base/cell-polygons.geojson.gz"

    # Check required files exist
    if [[ ! -f "$xenium_cells" ]]; then
        log_warn "  Comparison: cells.parquet not found in $region_dir — skipping."
        return 0
    fi
    if [[ ! -f "$proseg_cell_meta" ]]; then
        log_warn "  Comparison: proseg cell-metadata.csv.gz not found — skipping."
        return 0
    fi

    log_info "Comparing segmentations for $region_name ..."

    python3 - "$xenium_transcripts" "$xenium_cells" \
                "$proseg_cell_meta" "$proseg_tx_meta" \
                "$out_base" "$region_name" \
                "$xenium_boundaries" "$proseg_polygons" <<'PYEOF'
import sys, os, warnings, json, gzip
warnings.filterwarnings("ignore")

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("[WARN] pandas/numpy not available — skipping comparison.", file=sys.stderr)
    sys.exit(0)

try:
    import pyarrow.parquet as pq
    HAS_PARQUET = True
except ImportError:
    HAS_PARQUET = False

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection
    from matplotlib.lines import Line2D
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

(xenium_tx_path, xenium_cells_path,
 proseg_cell_path, proseg_tx_path,
 output_dir, region_name,
 xenium_boundaries_path, proseg_polygons_path) = sys.argv[1:9]

# ── Helper ────────────────────────────────────────────────────────
def read_parquet_safe(path):
    if not HAS_PARQUET:
        raise ImportError("pyarrow not installed")
    return pq.read_table(path).to_pandas()

def fmt(n):
    return f"{n:,.0f}"

def pct(num, denom):
    return f"{100 * num / denom:.1f}%" if denom > 0 else "N/A"

# ── Load original Xenium data ────────────────────────────────────
try:
    xen_cells = read_parquet_safe(xenium_cells_path)
except Exception as e:
    print(f"[WARN] Cannot read {xenium_cells_path}: {e}", file=sys.stderr)
    sys.exit(0)

xen_tx = None
if os.path.exists(xenium_tx_path):
    try:
        xen_tx = read_parquet_safe(xenium_tx_path)
    except Exception:
        pass

# ── Load proseg outputs ──────────────────────────────────────────
try:
    proseg_cells = pd.read_csv(proseg_cell_path)
except Exception as e:
    print(f"[WARN] Cannot read {proseg_cell_path}: {e}", file=sys.stderr)
    sys.exit(0)

proseg_tx = None
if os.path.exists(proseg_tx_path):
    try:
        proseg_tx = pd.read_csv(proseg_tx_path)
    except Exception:
        pass

# ── Compute metrics ──────────────────────────────────────────────
results = {"region": region_name}

# Cell counts
n_xen_cells = len(xen_cells)
n_proseg_cells = len(proseg_cells)
results["xenium_cells"] = n_xen_cells
results["proseg_cells"] = n_proseg_cells
results["cell_diff"] = n_proseg_cells - n_xen_cells
results["cell_diff_pct"] = 100 * (n_proseg_cells - n_xen_cells) / n_xen_cells if n_xen_cells > 0 else 0

# Transcripts per cell
xen_tx_col = None
for col in ["transcript_counts", "n_transcripts", "total_counts"]:
    if col in xen_cells.columns:
        xen_tx_col = col
        break

proseg_tx_col = None
for col in ["n_transcripts", "transcript_counts", "total_counts"]:
    if col in proseg_cells.columns:
        proseg_tx_col = col
        break

if xen_tx_col:
    xen_txpc = xen_cells[xen_tx_col].dropna()
    results["xenium_tx_per_cell_mean"] = float(xen_txpc.mean())
    results["xenium_tx_per_cell_median"] = float(xen_txpc.median())
    results["xenium_tx_per_cell_std"] = float(xen_txpc.std())

if proseg_tx_col:
    pro_txpc = proseg_cells[proseg_tx_col].dropna()
    results["proseg_tx_per_cell_mean"] = float(pro_txpc.mean())
    results["proseg_tx_per_cell_median"] = float(pro_txpc.median())
    results["proseg_tx_per_cell_std"] = float(pro_txpc.std())

# Cell area
xen_area_col = None
for col in ["cell_area", "area"]:
    if col in xen_cells.columns:
        xen_area_col = col
        break

proseg_area_col = None
for col in ["area", "cell_area"]:
    if col in proseg_cells.columns:
        proseg_area_col = col
        break

if xen_area_col:
    xen_area = xen_cells[xen_area_col].dropna()
    results["xenium_area_mean"] = float(xen_area.mean())
    results["xenium_area_median"] = float(xen_area.median())

if proseg_area_col:
    pro_area = proseg_cells[proseg_area_col].dropna()
    results["proseg_area_mean"] = float(pro_area.mean())
    results["proseg_area_median"] = float(pro_area.median())

# Transcript assignment rates
if xen_tx is not None:
    n_total_tx = len(xen_tx)
    cell_col = "cell_id" if "cell_id" in xen_tx.columns else None
    if cell_col:
        xen_assigned = xen_tx[cell_col].apply(
            lambda x: not (pd.isna(x) or str(x).upper() == "UNASSIGNED" or
                           (isinstance(x, (int, float)) and x < 0))
        ).sum()
        results["total_transcripts"] = n_total_tx
        results["xenium_assigned_tx"] = int(xen_assigned)
        results["xenium_assigned_pct"] = 100 * xen_assigned / n_total_tx if n_total_tx > 0 else 0

if proseg_tx is not None:
    n_proseg_total = len(proseg_tx)
    proseg_cell_col = None
    for col in ["cell", "cell_id", "assignment"]:
        if col in proseg_tx.columns:
            proseg_cell_col = col
            break
    if proseg_cell_col:
        proseg_assigned = proseg_tx[proseg_cell_col].apply(
            lambda x: not (pd.isna(x) or str(x).upper() == "UNASSIGNED" or
                           str(x).strip() == "" or
                           (isinstance(x, (int, float)) and x < 0))
        ).sum()
        results["proseg_assigned_tx"] = int(proseg_assigned)
        results["proseg_assigned_pct"] = 100 * proseg_assigned / n_proseg_total if n_proseg_total > 0 else 0

# ── Print report ─────────────────────────────────────────────────
sep = "─" * 60
print(f"\n{sep}")
print(f"  SEGMENTATION COMPARISON: {region_name}")
print(f"{sep}")

print(f"\n  {'Metric':<35} {'Xenium':>12} {'ProSeg':>12} {'Delta':>12}")
print(f"  {'─'*35} {'─'*12} {'─'*12} {'─'*12}")

print(f"  {'Cells':<35} {fmt(n_xen_cells):>12} {fmt(n_proseg_cells):>12}"
      f" {results['cell_diff']:>+12,.0f} ({results['cell_diff_pct']:+.1f}%)")

if xen_tx_col and proseg_tx_col:
    d = results["proseg_tx_per_cell_mean"] - results["xenium_tx_per_cell_mean"]
    print(f"  {'Transcripts/cell (mean)':<35}"
          f" {results['xenium_tx_per_cell_mean']:>12.1f}"
          f" {results['proseg_tx_per_cell_mean']:>12.1f}"
          f" {d:>+12.1f}")
    d = results["proseg_tx_per_cell_median"] - results["xenium_tx_per_cell_median"]
    print(f"  {'Transcripts/cell (median)':<35}"
          f" {results['xenium_tx_per_cell_median']:>12.1f}"
          f" {results['proseg_tx_per_cell_median']:>12.1f}"
          f" {d:>+12.1f}")

if xen_area_col and proseg_area_col:
    d = results["proseg_area_mean"] - results["xenium_area_mean"]
    print(f"  {'Cell area µm² (mean)':<35}"
          f" {results['xenium_area_mean']:>12.1f}"
          f" {results['proseg_area_mean']:>12.1f}"
          f" {d:>+12.1f}")
    d = results["proseg_area_median"] - results["xenium_area_median"]
    print(f"  {'Cell area µm² (median)':<35}"
          f" {results['xenium_area_median']:>12.1f}"
          f" {results['proseg_area_median']:>12.1f}"
          f" {d:>+12.1f}")

if "xenium_assigned_tx" in results:
    print(f"\n  {'Transcript Assignment':<35} {'Xenium':>12} {'ProSeg':>12}")
    print(f"  {'─'*35} {'─'*12} {'─'*12}")
    if "total_transcripts" in results:
        print(f"  {'Total transcripts':<35} {fmt(results['total_transcripts']):>12}")
    print(f"  {'Assigned transcripts':<35}"
          f" {fmt(results['xenium_assigned_tx']):>12}"
          f" {fmt(results.get('proseg_assigned_tx', 0)):>12}")
    print(f"  {'Assignment rate':<35}"
          f" {results['xenium_assigned_pct']:>11.1f}%"
          f" {results.get('proseg_assigned_pct', 0):>11.1f}%")

print(f"\n{sep}\n")

# ── Write TSV ────────────────────────────────────────────────────
tsv_path = os.path.join(output_dir, "segmentation_comparison.tsv")
pd.DataFrame([results]).to_csv(tsv_path, sep="\t", index=False)
print(f"  Comparison table saved: {tsv_path}")


# =====================================================================
# VISUAL COMPARISON: polygon borders side-by-side + overlay
# =====================================================================

if not HAS_MPL:
    print("[WARN] matplotlib not available — skipping visual comparison.")
    sys.exit(0)

if not HAS_PARQUET:
    print("[WARN] pyarrow not available — skipping visual comparison.")
    sys.exit(0)


# ── Parse Xenium cell_boundaries.parquet → dict{cell_id: [(x,y),...]} ─
def load_xenium_polygons(path):
    """
    Xenium cell_boundaries.parquet has columns:
      cell_id, vertex_x, vertex_y  (one row per vertex, grouped by cell_id)
    Returns dict mapping cell_id -> np.array of shape (N, 2).
    """
    df = pq.read_table(path, columns=["cell_id", "vertex_x", "vertex_y"]).to_pandas()
    polys = {}
    for cid, grp in df.groupby("cell_id"):
        polys[cid] = np.column_stack([grp["vertex_x"].values,
                                       grp["vertex_y"].values])
    return polys


# ── Parse proseg cell-polygons.geojson.gz → dict{cell_id: [ring, ...]} ─
def load_proseg_polygons(path):
    """
    Reads compressed GeoJSON FeatureCollection.
    Each Feature has geometry type Polygon or MultiPolygon.
    Returns dict mapping cell identifier -> list of np.array rings.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        gj = json.load(fh)

    polys = {}
    for feat in gj.get("features", []):
        props = feat.get("properties", {})
        # proseg may use "cell", "cell_id", or "id"
        cid = props.get("cell", props.get("cell_id", props.get("id", None)))
        geom = feat.get("geometry", {})
        gtype = geom.get("type", "")
        coords = geom.get("coordinates", [])

        rings = []
        if gtype == "Polygon":
            # coords = [exterior_ring, *holes]
            for ring in coords:
                rings.append(np.array(ring))
        elif gtype == "MultiPolygon":
            for poly in coords:
                for ring in poly:
                    rings.append(np.array(ring))
        if rings and cid is not None:
            polys[cid] = rings
    return polys


# ── Auto-select a zoom window around the densest centroid region ──
def find_dense_window(centroids_x, centroids_y, window_um=200, n_grid=50):
    """
    Tiles the tissue into a grid and picks the tile with the most
    centroids.  Returns (xmin, xmax, ymin, ymax) for the zoom window.
    """
    xmin_g, xmax_g = centroids_x.min(), centroids_x.max()
    ymin_g, ymax_g = centroids_y.min(), centroids_y.max()

    best_count = 0
    best_cx, best_cy = (xmin_g + xmax_g) / 2, (ymin_g + ymax_g) / 2

    step_x = (xmax_g - xmin_g) / n_grid
    step_y = (ymax_g - ymin_g) / n_grid
    half = window_um / 2

    for ix in range(n_grid):
        cx = xmin_g + (ix + 0.5) * step_x
        for iy in range(n_grid):
            cy = ymin_g + (iy + 0.5) * step_y
            mask = ((centroids_x >= cx - half) & (centroids_x <= cx + half) &
                    (centroids_y >= cy - half) & (centroids_y <= cy + half))
            c = mask.sum()
            if c > best_count:
                best_count = c
                best_cx, best_cy = cx, cy

    return (best_cx - half, best_cx + half,
            best_cy - half, best_cy + half)


# ── Filter polygons to a bounding box ────────────────────────────
def filter_xenium_polys(polys, bbox):
    """Keep only polygons whose centroid falls within bbox."""
    xmin, xmax, ymin, ymax = bbox
    out = {}
    for cid, verts in polys.items():
        cx, cy = verts[:, 0].mean(), verts[:, 1].mean()
        if xmin <= cx <= xmax and ymin <= cy <= ymax:
            out[cid] = verts
    return out


def filter_proseg_polys(polys, bbox):
    """Keep only polygons whose first ring centroid falls within bbox."""
    xmin, xmax, ymin, ymax = bbox
    out = {}
    for cid, rings in polys.items():
        r = rings[0]
        cx, cy = r[:, 0].mean(), r[:, 1].mean()
        if xmin <= cx <= xmax and ymin <= cy <= ymax:
            out[cid] = rings
    return out


# ── Draw polygon borders on an axis ─────────────────────────────
def draw_xenium_borders(ax, polys, bbox, color="#2166ac", lw=0.5, alpha=0.6):
    patches = []
    for cid, verts in polys.items():
        patches.append(MplPolygon(verts, closed=True))
    if patches:
        pc = PatchCollection(patches, facecolor="none",
                             edgecolor=color, linewidth=lw, alpha=alpha)
        ax.add_collection(pc)
    ax.set_xlim(bbox[0], bbox[1])
    ax.set_ylim(bbox[2], bbox[3])
    ax.set_aspect("equal")


def draw_proseg_borders(ax, polys, bbox, color="#b2182b", lw=0.5, alpha=0.6):
    patches = []
    for cid, rings in polys.items():
        for ring in rings:
            if len(ring) >= 3:
                patches.append(MplPolygon(ring[:, :2], closed=True))
    if patches:
        pc = PatchCollection(patches, facecolor="none",
                             edgecolor=color, linewidth=lw, alpha=alpha)
        ax.add_collection(pc)
    ax.set_xlim(bbox[0], bbox[1])
    ax.set_ylim(bbox[2], bbox[3])
    ax.set_aspect("equal")


# ── Gather centroid coordinates for zoom selection ───────────────
cx_col = None
for col in ["x_centroid", "centroid_x", "x", "X"]:
    if col in proseg_cells.columns:
        cx_col = col
        break
cy_col = None
for col in ["y_centroid", "centroid_y", "y", "Y"]:
    if col in proseg_cells.columns:
        cy_col = col
        break

if cx_col is None or cy_col is None:
    print("[WARN] Cannot determine centroid columns — skipping visual comparison.")
    sys.exit(0)

centroids_x = proseg_cells[cx_col].values.astype(float)
centroids_y = proseg_cells[cy_col].values.astype(float)

# ── Check polygon files exist ────────────────────────────────────
if not os.path.exists(xenium_boundaries_path):
    print(f"[WARN] {xenium_boundaries_path} not found — skipping visual comparison.")
    sys.exit(0)
if not os.path.exists(proseg_polygons_path):
    print(f"[WARN] {proseg_polygons_path} not found — skipping visual comparison.")
    sys.exit(0)

print(f"  Loading polygon data for visual comparison...")

try:
    xen_polys = load_xenium_polygons(xenium_boundaries_path)
except Exception as e:
    print(f"[WARN] Failed to load Xenium boundaries: {e}")
    sys.exit(0)

try:
    proseg_polys = load_proseg_polygons(proseg_polygons_path)
except Exception as e:
    print(f"[WARN] Failed to load proseg polygons: {e}")
    sys.exit(0)

print(f"  Xenium polygons loaded: {len(xen_polys):,}")
print(f"  ProSeg polygons loaded: {len(proseg_polys):,}")

# ── Compute zoom windows ────────────────────────────────────────
#    Two zoom levels: 200 µm (cellular detail) and 500 µm (neighbourhood)
zoom_specs = [
    ("detail_200um", 200),
    ("neighbourhood_500um", 500),
]

for zoom_tag, window_um in zoom_specs:
    bbox = find_dense_window(centroids_x, centroids_y, window_um=window_um)
    xen_z = filter_xenium_polys(xen_polys, bbox)
    pro_z = filter_proseg_polys(proseg_polys, bbox)

    if len(xen_z) < 3 and len(pro_z) < 3:
        print(f"  [{zoom_tag}] Too few cells in window — skipping.")
        continue

    print(f"  [{zoom_tag}] Xenium: {len(xen_z)} cells | ProSeg: {len(pro_z)} cells")

    # 3-panel figure: Xenium | ProSeg | Overlay
    fig, axes = plt.subplots(1, 3, figsize=(21, 7), dpi=150)

    xen_color = "#2166ac"
    pro_color = "#b2182b"

    # Panel 1: Xenium original
    draw_xenium_borders(axes[0], xen_z, bbox, color=xen_color, lw=0.6)
    axes[0].set_title(f"Xenium Original ({len(xen_z)} cells)",
                      fontsize=12, fontweight="bold", color=xen_color)
    axes[0].set_xlabel("x (µm)", fontsize=9)
    axes[0].set_ylabel("y (µm)", fontsize=9)
    axes[0].tick_params(labelsize=7)

    # Panel 2: ProSeg
    draw_proseg_borders(axes[1], pro_z, bbox, color=pro_color, lw=0.6)
    axes[1].set_title(f"ProSeg ({len(pro_z)} cells)",
                      fontsize=12, fontweight="bold", color=pro_color)
    axes[1].set_xlabel("x (µm)", fontsize=9)
    axes[1].set_ylabel("y (µm)", fontsize=9)
    axes[1].tick_params(labelsize=7)

    # Panel 3: Overlay
    draw_xenium_borders(axes[2], xen_z, bbox, color=xen_color, lw=0.5, alpha=0.5)
    draw_proseg_borders(axes[2], pro_z, bbox, color=pro_color, lw=0.5, alpha=0.5)
    legend_handles = [
        Line2D([0], [0], color=xen_color, lw=1.5, label="Xenium"),
        Line2D([0], [0], color=pro_color, lw=1.5, label="ProSeg"),
    ]
    axes[2].legend(handles=legend_handles, loc="upper right", fontsize=9)
    axes[2].set_title("Overlay", fontsize=12, fontweight="bold")
    axes[2].set_xlabel("x (µm)", fontsize=9)
    axes[2].set_ylabel("y (µm)", fontsize=9)
    axes[2].tick_params(labelsize=7)

    fig.suptitle(f"{region_name} — Segmentation Comparison ({window_um} µm window)",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()

    fig_path = os.path.join(output_dir, f"segmentation_comparison_{zoom_tag}.pdf")
    fig.savefig(fig_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Visual comparison saved: {fig_path}")


# ── Distribution comparison figure ───────────────────────────────
#    Side-by-side histograms: transcripts/cell + cell area
n_dist_panels = 0
if xen_tx_col and proseg_tx_col:
    n_dist_panels += 1
if xen_area_col and proseg_area_col:
    n_dist_panels += 1

if n_dist_panels > 0:
    fig, axes = plt.subplots(1, n_dist_panels, figsize=(7 * n_dist_panels, 5), dpi=150)
    if n_dist_panels == 1:
        axes = [axes]
    pi = 0

    if xen_tx_col and proseg_tx_col:
        ax = axes[pi]; pi += 1
        bins = np.linspace(0,
                           max(xen_txpc.quantile(0.99), pro_txpc.quantile(0.99)),
                           60)
        ax.hist(xen_txpc, bins=bins, alpha=0.55, color="#2166ac",
                label=f"Xenium (med={xen_txpc.median():.0f})", density=True)
        ax.hist(pro_txpc, bins=bins, alpha=0.55, color="#b2182b",
                label=f"ProSeg (med={pro_txpc.median():.0f})", density=True)
        ax.set_xlabel("Transcripts per cell", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.set_title("Transcripts per Cell", fontsize=12, fontweight="bold")
        ax.legend(fontsize=9)

    if xen_area_col and proseg_area_col:
        ax = axes[pi]; pi += 1
        bins = np.linspace(0,
                           max(xen_area.quantile(0.99), pro_area.quantile(0.99)),
                           60)
        ax.hist(xen_area, bins=bins, alpha=0.55, color="#2166ac",
                label=f"Xenium (med={xen_area.median():.0f})", density=True)
        ax.hist(pro_area, bins=bins, alpha=0.55, color="#b2182b",
                label=f"ProSeg (med={pro_area.median():.0f})", density=True)
        ax.set_xlabel("Cell area (µm²)", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.set_title("Cell Area Distribution", fontsize=12, fontweight="bold")
        ax.legend(fontsize=9)

    fig.suptitle(f"{region_name} — Distribution Comparison",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()

    dist_path = os.path.join(output_dir, "segmentation_comparison_distributions.pdf")
    fig.savefig(dist_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  Distribution plot saved: {dist_path}")

print(f"\n  All comparisons complete for {region_name}.")

PYEOF

    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        log_warn "  Comparison script exited with code $exit_code"
    fi
    return 0
}

# ============================================================================
# ============================================================================
# Build a dense pixel-level mask from Xenium's cell_boundaries.parquet for use
# as proseg's --cellpose-masks prior.
#
# Arguments:
#   $1 = region directory containing cell_boundaries.parquet + transcripts.parquet
#   $2 = output directory (where to write mask + params)
#   $3 = mask scale (microns per pixel, e.g. 0.5)
#
# On success, writes:
#   <out>/xenium_polygons_mask.npy.gz      (uint32 mask, gzipped)
#   <out>/cellpose_params.env              (CELLPOSE_X_MIN, CELLPOSE_Y_MIN)
#
# The mask values are 1..N (0 = background) where the integer corresponds to
# label_id in Xenium's cell_boundaries.parquet, matching the row order of
# cells.parquet. Mapping back to hex cell_ids is preserved in the proseg output.
# ============================================================================
build_xenium_mask() {
    local region_dir="$1"
    local out_dir="$2"
    local scale="$3"

    log_info "  Rasterising Xenium polygons → mask at ${scale} µm/px"

    python3 - "$region_dir" "$out_dir" "$scale" <<'PYEOF'
import sys
import gzip
import numpy as np
import pandas as pd
from pathlib import Path
from shapely.geometry import Polygon
from rasterio.features import rasterize
from rasterio.transform import Affine

region_dir = Path(sys.argv[1])
out_dir    = Path(sys.argv[2])
scale      = float(sys.argv[3])

bounds_path = region_dir / "cell_boundaries.parquet"
tx_path     = region_dir / "transcripts.parquet"
if not bounds_path.exists():
    sys.exit(f"ERROR: {bounds_path} not found")
if not tx_path.exists():
    sys.exit(f"ERROR: {tx_path} not found")

# Spatial extent comes from transcripts (the coordinate frame proseg operates in)
tx = pd.read_parquet(tx_path, columns=["x_location", "y_location"])
buf = 5.0  # micron buffer so cells at the edge aren't clipped
x_min = float(tx["x_location"].min()) - buf
x_max = float(tx["x_location"].max()) + buf
y_min = float(tx["y_location"].min()) - buf
y_max = float(tx["y_location"].max()) + buf

# proseg's clap-rs argument parser rejects values starting with '-<digit>' (it
# interprets them as short flags), so we cannot pass negative origins via
# --cellpose-{x,y}-transform. Clamp the mask origin to be non-negative.
# Polygons with vertices at slightly negative coordinates will have those edges
# clipped during rasterisation, but the bulk of each cell remains in the mask.
_x_clamped = max(0.0, x_min)
_y_clamped = max(0.0, y_min)
if _x_clamped != x_min or _y_clamped != y_min:
    print(f"  Origin clamped to ≥0 to avoid proseg parser bug: "
          f"({x_min:.3f}, {y_min:.3f}) → ({_x_clamped:.3f}, {_y_clamped:.3f})",
          flush=True)
x_min, y_min = _x_clamped, _y_clamped

width  = int(np.ceil((x_max - x_min) / scale))
height = int(np.ceil((y_max - y_min) / scale))
mb = (width * height * 4) / 1e6
print(f"  Extent: x∈[{x_min:.1f},{x_max:.1f}] y∈[{y_min:.1f},{y_max:.1f}] µm",
      flush=True)
print(f"  Mask:   {height}×{width} px → {mb:,.1f} MB uncompressed", flush=True)

# Load polygon vertices; use label_id directly (it's already 1..N uint32-ready)
bounds = pd.read_parquet(bounds_path)
required = {"cell_id", "vertex_x", "vertex_y", "label_id"}
missing = required - set(bounds.columns)
if missing:
    sys.exit(f"ERROR: cell_boundaries.parquet missing columns: {missing}")

shapes = []
invalid = 0
for (cid, lab), grp in bounds.groupby(["cell_id", "label_id"], sort=False):
    pts = list(zip(grp["vertex_x"].values, grp["vertex_y"].values))
    if len(pts) < 3:
        invalid += 1; continue
    p = Polygon(pts)
    if not p.is_valid:
        p = p.buffer(0)  # repair self-intersection
    if not p.is_valid or p.is_empty:
        invalid += 1; continue
    shapes.append((p, int(lab)))

print(f"  Polygons: {len(shapes):,} valid, {invalid} skipped", flush=True)
if not shapes:
    sys.exit("ERROR: no valid polygons to rasterise")

# Affine maps (col, row) → (x_world, y_world)
# pixel (0,0) sits at world (x_min, y_min); both axes grow positively with index
transform = Affine(scale, 0, x_min, 0, scale, y_min)

mask = rasterize(
    shapes,
    out_shape=(height, width),
    transform=transform,
    dtype=np.uint32,
    fill=0,
)

occupied = (mask > 0).sum() / mask.size * 100
n_cells  = len(np.unique(mask)) - 1  # subtract background
print(f"  Mask occupied: {occupied:.1f}% — {n_cells:,} cell labels visible",
      flush=True)

out_dir.mkdir(parents=True, exist_ok=True)
mask_path = out_dir / "xenium_polygons_mask.npy.gz"
with gzip.open(mask_path, "wb") as f:
    np.save(f, mask)
sz_mb = mask_path.stat().st_size / 1e6
print(f"  Wrote {mask_path.name} ({sz_mb:,.1f} MB compressed)", flush=True)

# Emit shell-sourcable params for the proseg invocation
env_path = out_dir / "cellpose_params.env"
env_path.write_text(
    f"CELLPOSE_X_MIN={x_min}\n"
    f"CELLPOSE_Y_MIN={y_min}\n"
    f"CELLPOSE_WIDTH={width}\n"
    f"CELLPOSE_HEIGHT={height}\n"
    f"CELLPOSE_N_CELLS={n_cells}\n"
)
print(f"  Wrote {env_path.name}", flush=True)
PYEOF

    local rc=$?
    if [[ $rc -ne 0 ]]; then
        log_error "Mask rasterisation failed for $region_dir (exit $rc)"
        return 1
    fi
    return 0
}

# ============================================================================
# Process a single region
# ============================================================================
process_region() {
    local parquet_path="$1"
    local region_dir
    local region_name
    local out_base

    region_dir="$(dirname "$parquet_path")"
    region_name="$(basename "$region_dir")"

    local safe_name
    safe_name="$(echo "$region_name" | tr ' ' '_')"

    out_base="${OUTPUT_DIR}/${safe_name}"

    log_info "========================================"
    log_info "Processing: $region_name"
    log_info "  Input:  $parquet_path"
    log_info "  Output: $out_base"
    log_info "========================================"

    if $DRY_RUN; then
        log_info "[DRY RUN] Would process $region_name"
        return 0
    fi

    mkdir -p "$out_base"
    mkdir -p "$out_base/xenium_explorer"

    # ------------------------------------------------------------------
    # Step 1: Run proseg
    # ------------------------------------------------------------------
    local proseg_cmd=(
        "$PROSEG_BIN"
        --xenium
        --overwrite

        # ── Always-on outputs ─────────────────────────────────────────
        --output-spatialdata         "$out_base/proseg-output.zarr"
        --output-counts              "$out_base/counts.mtx.gz"
        --output-cell-metadata       "$out_base/cell-metadata.csv.gz"
        --output-gene-metadata       "$out_base/gene-metadata.csv.gz"
        --output-transcript-metadata "$out_base/transcript-metadata.csv.gz"
        --output-cell-polygons       "$out_base/cell-polygons.geojson.gz"
    )

    # ── Optional output: expected (non-integer) count matrix ──────────
    if $OUTPUT_EXPECTED_COUNTS; then
        proseg_cmd+=(
            --output-expected-counts "$out_base/expected-counts.mtx.gz"
        )
    fi

    # ── Optional output: cell-by-gene Poisson rate parameters ─────────
    if $OUTPUT_RATES; then
        proseg_cmd+=(
            --output-rates "$out_base/rates.csv.gz"
        )
    fi

    # ── Optional output: per-z-layer non-overlapping polygons ─────────
    if $OUTPUT_POLYGON_LAYERS; then
        proseg_cmd+=(
            --output-cell-polygon-layers "$out_base/cell-polygons-layers.geojson.gz"
        )
    fi

    # ── Optional: omit transcript positions from zarr ─────────────────
    if $EXCLUDE_TRANSCRIPTS; then
        proseg_cmd+=( --exclude-spatialdata-transcripts )
    fi

    # ── Threads ───────────────────────────────────────────────────────
    if [[ -n "$THREADS" ]]; then
        proseg_cmd+=(--nthreads "$THREADS")
    fi

    # ── All model / sampling args forwarded verbatim ──────────────────
    if [[ ${#EXTRA_PROSEG_ARGS[@]} -gt 0 ]]; then
        proseg_cmd+=("${EXTRA_PROSEG_ARGS[@]}")
    fi

    # ── Xenium-mask prior (rasterise polygons → cellpose-format mask) ─
    if $USE_XENIUM_MASK; then
        if ! build_xenium_mask "$region_dir" "$out_base" "$MASK_SCALE"; then
            log_error "  Skipping $region_name due to mask build failure"
            return 1
        fi
        # Source CELLPOSE_X_MIN / CELLPOSE_Y_MIN written by the Python helper
        # shellcheck source=/dev/null
        source "$out_base/cellpose_params.env"
        # Proseg requires either --cellpose-scale (uniform scale only)
        # OR both --cellpose-{x,y}-transform (general affine); they are
        # mutually exclusive. We use the transforms because our mask
        # origin is not at (0, 0) µm — the scale is encoded as the 'a'
        # coefficient of each transform.
        proseg_cmd+=(
            --cellpose-masks       "$out_base/xenium_polygons_mask.npy.gz"
            --cellpose-x-transform "$MASK_SCALE" 0 "$CELLPOSE_X_MIN"
            --cellpose-y-transform 0 "$MASK_SCALE" "$CELLPOSE_Y_MIN"
        )
        log_info "  Using Xenium mask prior: ${CELLPOSE_N_CELLS} cells, "\
"origin (${CELLPOSE_X_MIN}, ${CELLPOSE_Y_MIN}) µm"
    fi

    # ── Input file (must be last) ──────────────────────────────────────
    proseg_cmd+=("$parquet_path")

    log_info "Running proseg..."
    log_info "  Command: ${proseg_cmd[*]}"

    local start_time
    start_time=$(date +%s)

    if "${proseg_cmd[@]}" 2>&1 | tee "$out_base/proseg.log"; then
        local end_time elapsed
        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))
        log_ok "proseg completed for $region_name in ${elapsed}s"
    else
        local end_time elapsed
        end_time=$(date +%s)
        elapsed=$(( end_time - start_time ))
        log_error "proseg FAILED for $region_name after ${elapsed}s"
        log_error "  Check log: $out_base/proseg.log"
        return 1
    fi

    # ------------------------------------------------------------------
    # Step 2: Convert for Xenium Explorer (proseg-to-baysor)
    # ------------------------------------------------------------------
    if [[ -n "$PROSEG_TO_BAYSOR" ]]; then
        log_info "Converting to Xenium Explorer format..."

        local baysor_cmd=(
            "$PROSEG_TO_BAYSOR"
            "$out_base/proseg-output.zarr"
            --output-transcript-metadata "$out_base/xenium_explorer/transcript-metadata.csv"
            --output-cell-polygons       "$out_base/xenium_explorer/cell-polygons.geojson"
        )

        if "${baysor_cmd[@]}" 2>&1 | tee -a "$out_base/proseg.log"; then
            log_ok "Xenium Explorer conversion complete for $region_name"

            cat > "$out_base/xenium_explorer/import_to_xenium_ranger.sh" <<XENIUM_EOF
#!/usr/bin/env bash
# Import proseg segmentation into Xenium Explorer via xeniumranger.
# Edit --xenium-bundle to point to the ORIGINAL Xenium output bundle.

xeniumranger import-segmentation \\
    --id="${safe_name}_proseg" \\
    --xenium-bundle="${region_dir}" \\
    --transcript-assignment="$out_base/xenium_explorer/transcript-metadata.csv" \\
    --viz-polygons="$out_base/xenium_explorer/cell-polygons.geojson" \\
    --units=microns \\
    --localcores=32 \\
    --localmem=128
XENIUM_EOF
            chmod +x "$out_base/xenium_explorer/import_to_xenium_ranger.sh"
            log_info "  Xenium Ranger import script: $out_base/xenium_explorer/import_to_xenium_ranger.sh"
        else
            log_warn "proseg-to-baysor conversion failed for $region_name"
        fi
    fi

    # ------------------------------------------------------------------
    # Step 3: Compare original vs proseg segmentation
    # ------------------------------------------------------------------
    if $COMPARE_SEG; then
        compare_segmentation "$region_dir" "$out_base" "$region_name"
    fi

    log_ok "Region complete: $region_name"
    return 0
}

# ============================================================================
# Execute across all regions
# ============================================================================
TOTAL=${#PARQUET_FILES[@]}
SUCCESS=0
FAILED=0
FAILED_REGIONS=()

log_info "Starting proseg pipeline"
log_info "  Mode:         $([ "$PARALLEL_JOBS" -gt 1 ] && echo "parallel ($PARALLEL_JOBS jobs)" || echo "sequential")"
log_info "  Regions:      $TOTAL"
log_info "  Output root:  $OUTPUT_DIR/"
if $USE_XENIUM_MASK; then
    log_info "  Prior:        Xenium polygon mask (--cellpose-masks, ${MASK_SCALE} µm/px)"
else
    log_info "  Prior:        transcript-table cell_id column (default)"
fi
if [[ ${#EXTRA_PROSEG_ARGS[@]} -gt 0 ]]; then
    log_info "  Extra args:   ${EXTRA_PROSEG_ARGS[*]}"
fi
echo ""

if [[ "$PARALLEL_JOBS" -le 1 ]]; then
    for i in "${!PARQUET_FILES[@]}"; do
        log_info "Region $(( i + 1 )) / $TOTAL"
        if process_region "${PARQUET_FILES[$i]}"; then
            (( SUCCESS++ )) || true
        else
            (( FAILED++ )) || true
            FAILED_REGIONS+=("$(basename "$(dirname "${PARQUET_FILES[$i]}")")")
        fi
        echo ""
    done
else
    export INPUT_DIR OUTPUT_DIR THREADS DRY_RUN PROSEG_BIN PROSEG_TO_BAYSOR
    export OUTPUT_EXPECTED_COUNTS OUTPUT_RATES OUTPUT_POLYGON_LAYERS EXCLUDE_TRANSCRIPTS
    export COMPARE_SEG
    export USE_XENIUM_MASK MASK_SCALE
    export -a EXTRA_PROSEG_ARGS
    export -f process_region compare_segmentation log_info log_ok log_warn log_error build_xenium_mask

    if command -v parallel &>/dev/null; then
        log_info "Using GNU parallel with $PARALLEL_JOBS jobs"
        printf '%s\n' "${PARQUET_FILES[@]}" | \
            parallel -j "$PARALLEL_JOBS" --halt soon,fail=1 \
            process_region {}
    else
        log_warn "GNU parallel not found. Falling back to background jobs."
        log_warn "  Install with: sudo apt install parallel"

        local_running=0
        for parquet_path in "${PARQUET_FILES[@]}"; do
            process_region "$parquet_path" &
            (( local_running++ )) || true
            if (( local_running >= PARALLEL_JOBS )); then
                wait -n
                (( local_running-- )) || true
            fi
        done
        wait
    fi

    # Tally results for parallel mode
    for f in "${PARQUET_FILES[@]}"; do
        region_dir="$(dirname "$f")"
        region_name="$(basename "$region_dir")"
        safe_name="$(echo "$region_name" | tr ' ' '_')"
        if [[ -d "${OUTPUT_DIR}/${safe_name}/proseg-output.zarr" ]]; then
            (( SUCCESS++ )) || true
        else
            (( FAILED++ )) || true
            FAILED_REGIONS+=("$region_name")
        fi
    done
fi

# ============================================================================
# Summary
# ============================================================================
echo ""
log_info "========================================"
log_info "Pipeline complete"
log_info "  Total:      $TOTAL"
log_ok   "  Succeeded:  $SUCCESS"
if [[ $FAILED -gt 0 ]]; then
    log_error "  Failed:     $FAILED"
    for r in "${FAILED_REGIONS[@]}"; do
        log_error "    - $r"
    done
fi
log_info "  Output at:  $OUTPUT_DIR/"
log_info "========================================"
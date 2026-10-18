# Passive thermography corrosion assessment (MATLAB)

Passive thermography for detecting and sizing corrosion on steel plates, as a basis for wind turbine tower inspection. No heat source is used: the plate is imaged with a thermal IR camera as it heats and cools naturally, and corroded areas show up as uneven temperature across the surface.

A test plate is compared against a bare-steel reference plate at the same plate temperature (2 degC bands), using the spread (std) of the plate temperature. The result is a verdict (bare-like, uniform or patchy), an estimate of the affected area, and a class for each band.

## Pipeline order

```
IRSoft export  ->  check_roi_fit  ->  corrosion_assess
(IV_#####.xlsx)    (per folder)       (test vs reference)
```

### 0. Data
One folder per plate under `root`, each holding thermal frames exported from IRSoft as `IV_#####.xlsx` (480x640, degC).
- Reference folder: bare steel (e.g. `BARE NEW`)
- Test folder: the plate being assessed (e.g. `MIXEDB`)

### 1. `check_roi_fit`
Run on the reference folder and on every test folder.
- Cuts the plate region out of each frame and checks the edges for background or cold pixels.
- Flags bad frames as `CHECK: ...`, good frames as `OK`.
- Assigns each frame to a 2 degC plate-temperature band.
- Writes `stats\frame_bands_<foldername>.csv`.

### 2. `corrosion_assess`
Run once per test folder, after step 1.
- Reads the `frame_bands_<foldername>.csv` files and drops every frame flagged `CHECK` (test and reference).
- Matches test and reference frames of the same band, then gives per band:
  - std-ratio verdict: BARE-LIKE / UNIFORM / PATCHY
  - hot-area size (steel-referenced)
  - class: LOCALISED / PARTIAL / WIDESPREAD / UNIFORM / UNEVEN / BARE-LIKE / MIXED
- Output in `stats\`: `assess_<testName>.csv`, `_frames.csv`, `_run.txt`, `.png` (std ratio vs temperature plot).

## Helper functions (called by the two scripts above)
- `extract_plate_roi`: plate ROI cut and masking
- `load_thermal_matrix`: reads and validates one frame
- `compute_indicators_multi`: std ratio of test vs reference frames

## Setup
Set `root` at the top of each script (two lines, comment out the one you are not using). All plate folders and `stats\` sit under it.

## Additional scripts 
Olptional scripts, not needed to assess a new plate: `batch_compare_bands_*`, `background_check`, `thermal_anomaly_map`, `compute_indicators`, `uniformity_index`, and the photo-validation scripts (`tune_rust`, `refine_rust`, `rust_colour_rel`, `validate_with_photo`).
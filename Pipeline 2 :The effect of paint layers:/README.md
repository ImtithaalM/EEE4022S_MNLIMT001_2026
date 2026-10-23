# Passive thermography: Effect of paint layers 

Test 2 of the passive thermography corrosion work. Rusted steel plates were coated with one or two layers of paint, and the question is whether the corrosion can still be detected and sized through the paint. The method is the same as Pipeline 1 (corrosion_assess), but each coated plate is compared with a clean coated reference that has the same number of coats, because paint changes the apparent temperature of the whole surface.

The result per plate and 2 degC band is a verdict (bare-like, uniform or patchy), an estimate of the corroded area, and a class, plus a roll-up across the good bands for each plate.

## Pipeline 2 order

```
IRSoft export  ->  check_roi_fit  ->  painted_assess
(IV_#####.xlsx)    (per folder)       (all coated plates, each vs its own reference)
```

### 0. Data
One folder per plate under `root`, each holding thermal frames exported from IRSoft as `IV_#####.xlsx` (480x640, degC).
- Reference folders: clean coated steel, `BARE PAINT` (1 coat) and `BARE PAINT2` (2 coats)
- Test folders (rusted, then coated):
  - 1 coat: `PAINT15`, `PAINT50`, `PAINT100`
  - 2 coats: `PAINT152`, `PAINT502`, `PAINT1002`
  - the number is the approximate rust share (light, half, full) and a trailing 2 means two coats

### 1. `check_roi_fit`
Run on both reference folders and on every test folder.
- Cuts the plate region out of each frame and checks the edges for background or cold pixels.
- Flags bad frames as `CHECK: ...`, good frames as `OK`.
- Assigns each frame to a 2 degC plate-temperature band.
- Writes `stats\frame_bands_<foldername>.csv`.

### 2. `painted_assess`
Run once. It loops over the plates listed in `plates` (name, reference folder, photo rust %).
- Reads the `frame_bands_<foldername>.csv` files and drops every frame flagged `CHECK` (test and reference).
- Matches test and reference frames of the same band (same coat count), then gives per band:
  - std-ratio verdict: BARE-LIKE / UNIFORM / PATCHY
  - hot-area size (steel-referenced)
  - class: LOCALISED / PARTIAL / WIDESPREAD / UNIFORM / UNEVEN / BARE-LIKE / MIXED
- Roll-up per plate over the good bands only (at least 3 reference and 3 test frames, and session closeness of at most 3 degC): median and range of the std ratio and area, most common class, and the photo rust %.
- Output in `stats\`: `assess_<plate>.csv`, `assess_<plate>_frames.csv`, `assess_<plate>_run.txt`, and `painted_assess_summary.csv`.

## Helper functions (called by the two scripts above)
- `extract_plate_roi`: plate ROI cut and masking
- `load_thermal_matrix`: reads and validates one frame
- `compute_indicators_multi`: std ratio of test vs reference frames

## Setup
Set `root` at the top of `painted_assess` (two lines, comment out the one you are not using). All plate folders and `stats\` sit under it. To add a plate, add one entry to `plates` with its folder, its reference and the photo rust %.

## Settings worth knowing
Same detection and sizing settings as Pipeline 1 (2 degC bands, margin 0.05, area cut-offs for localised and widespread). The roll-up thresholds `goodMinFrames` and `goodCloseness` decide which bands count as good.
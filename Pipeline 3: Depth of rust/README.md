# Passive thermography: Depth of rust 

Test 3 of the passive thermography corrosion work. Two steel plates (2 mm and 4 mm) carry rusted spots of different thickness. The micrometer thickness of each spot (with rust minus after the rust was removed) is compared with the temperature difference dT between the spot and the rest of the plate, in each temperature band. The question is whether thicker rust gives a larger dT.

## Pipeline 3 order

```
IRSoft export  ->  check_depth_roi_fit  ->  analyse_depth
(IV_#####.xlsx)    (per plate folder,       (dT vs rust thickness
                    uses depth_geometry)     per band)
```

### 0. Data
- One folder per plate under `baseDir`: `2mm new` and `4mm new`, each holding thermal frames exported from IRSoft as `IV_#####.xlsx` (480x640, degC).
- Micrometer readings, three per spot, before (with rust) and after (bare steel), typed into `analyse_depth` under Settings.
- Active spots: 2 mm plate TM and BM, 4 mm plate TR2 and B. The original four spots on each plate went shiny and are not analysed (their readings are kept as comments).

### 1. `depth_geometry`
Called by `check_depth_roi_fit`, not run on its own. `g = depth_geometry('2mm')` or `('4mm')` returns:
- the fixed plate outline (four lines, hard-coded from IRSoft line profiles; redo only if the plate or camera moves)
- the depth rectangles, which of them are active, and the region kept out of the reference
- the inset, hot-speck mask and reference-region settings

### 2. `check_depth_roi_fit`
Run on each plate folder.
- Cuts the plate region out of each frame with the fixed outline, masks hot specks, and checks the fit.
- Flags bad frames, and assigns each good frame to a 4 degC band (22-26, 28-32, 38-42, 48-52) or `OUT`.
- For each active rectangle, writes dT: the rectangle median minus the plate reference median, where the reference is the plate outside all rectangles and the excluded region.
- Writes `stats\depth_frame_bands_<foldername>.csv` with `status`, `band` and one `<rect>_dT` column per active rectangle.

### 3. `analyse_depth`
Run after step 2 for both plates.
- Rust thickness per spot = mean of before readings minus mean of after readings, with a standard error from the six readings.
- Keeps frames with status `OK` and a band that is not `OUT`.
- Per band and spot: mean dT and its standard deviation over the usable frames.
- Two figures: dT against rust thickness with one line per band (both plates), and one panel per band.
- Output in `stats\`: `depth_vs_thickness_new.csv`. `analyse_depth_combined` adds the original four spots.

## Helper functions
- `load_thermal_matrix`: reads and validates one frame
- `depth_geometry`: plate outline and rectangles (above)

## Setup
Set `baseDir` at the top of `analyse_depth` (two lines, comment out the one you are not using). The plate folders and `stats\` sit under it. If you re-measure a spot, change its before and after rows in `plates`. If the camera or plate moves, redo the line profiles in `depth_geometry` first.

## Notes
- The micrometer anvils rest on the high points of a rough rusted surface, so the measured rust thickness is likely underestimated.
- The result is a trend in dT with thickness, not a calibrated depth measurement.
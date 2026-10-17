# EEE4022S_MNLIMT001_2026

Code, scripts and data-processing pipelines for my EEE4022S final-year capstone at the University of Cape Town (BSc Eng, Mechatronics / Electrical Engineering).

**Thesis title:** Passive thermography as corrosion detection for wind turbine tower

**Student:** Imtithaal Manuel (MNLIMT001)
**Supervisor:** Prof. Khan

---

## Project overview

Wind turbine towers are steel structures that corrode over time, and inspecting them is slow, costly and dangerous. This project investigates whether **passive thermography** (using the natural temperature variation of the tower, e.g. solar heating, instead of an active heat source) can detect corrosion on steel from thermal camera images alone.

The project began as drone-based inspection (AK-03, image-processing only). After evaluating a range of NDT techniques against the drone/no-contact constraint, it settled on passive thermography as the primary method. Testing is done in a controlled lab on mild steel plates (~100 x 75 x 10 mm, plus 2 mm and 4 mm plates). A heat gun is used only to set the plate temperature as a stand-in for field solar heating.

### Equipment

- Testo 882 thermal camera (320 x 240, NETD < 50 mK)
- Sentry ST642 non-contact infrared thermometer

### Test variables

- Steel thickness (2 / 4 / 10 mm)
- Corrosion level (none / light / 50% / 100%)
- Temperature (20-50 °C)
- Coatings (1-coat and 2-coat paint)
- Depth of rust

---

## Repository structure

```
EEE4022S_MNLIMT001_2026/
├── README.md
├── pipeline1_corrosion_detection/   # bare vs corroded (light / HALF / FULL) comparison
├── pipeline2_coated_plates/         # painted plates vs painted bare reference
├── pipeline3_depth_of_rust/         # corrosion thickness estimation
├── data_checklist/                  # data checklist / logs of captured datasets
└── docs/                            # experimental setup procedure, notes
```

> Raw thermal images/data are not stored in this repo (size). See [Data](#data).

---

## The three pipelines

### Pipeline 1: Corrosion detection (bare vs corroded)

**Question:** Can thermal images tell a corroded plate apart from a bare one?

- Compares a bare reference plate against corroded plates (light, 50%, FULL).
- Extracts a region of interest (ROI) per plate, then computes a standard-deviation ratio between corroded and bare to give a detection verdict.
- `rev2_demo.m` is the quick bare vs FULL vs light std-ratio verdict demo.
- Note: the FULL dataset was retaken unpainted. An earlier painted version looked artificially uniform and is treated as invalid for this comparison.

**Inputs:** thermal images/exports from the Testo 882 per plate and temperature
**Outputs:** per-plate statistics, std-ratio verdicts, comparison plots

### Pipeline 2: Coated plates (paint)

**Question:** Does a paint layer hide corrosion from passive thermography?

- Same approach as Pipeline 1, but compares painted corroded plates against a painted bare reference.
- Folder naming:

| Folder | Plate | Coats |
|---|---|---|
| `PAINT50` | HALF plate | 1 |
| `PAINT502` | HALF plate | 2 |
| `PAINT100` | FULL plate | 1 |
| `PAINT1002` | FULL plate | 2 |
| `BARE PAINT` | Painted bare reference | n/a |

- The HALF plate is physically smaller than the others, so each plate needs its own ROI outline.
- Currently running the 1-coat plates (PAINT50, PAINT100) first. The 2-coat plates are parked until their ROI outlines exist.
- Per supervisor guidance, the 1-coat vs 2-coat trend is reported as a result/discussion point, with no further coating work.

### Pipeline 3: Depth of rust

**Question:** Can thermal response indicate how thick the corrosion layer is?

- Uses fully corroded 2 mm and 4 mm plates, with corrosion removed at several spots per plate to create different depths.
- Ground truth: micrometer thickness measurements, with each spot stripped to bare steel to get the corrosion thickness.
- Only the two newer spots per plate are analysed, because the original spots went shiny and read low.
- Built as a separate MATLAB pipeline: 2 mm first, then 4 mm.


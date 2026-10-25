# Data

## `sint800.tsv`: the Sint-800 record (Guyodo & Valet 1999)

The composite record of relative palaeointensity of the past 800 kyr and its calibration to
the virtual axial dipole moment (VADM), from Y. Guyodo & J.-P. Valet, "Global changes in
intensity of the Earth's magnetic field during the past 800 kyr", Nature 399, 249–252
(1999), doi:10.1038/20420. This is the record analysed by Hoyng, Schmitt & Ossendrijver
(2002) (Tier 6 of `dev/VALIDATION.md`; `src/sint800.jl`).

Provenance: the third Supplementary Information file of the article,
`41586_1999_BF20420_MOESM3_ESM.rtf` (retrieved 2026-10-24 from
https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2F20420/MediaObjects/41586_1999_BF20420_MOESM3_ESM.rtf,
SHA-256 `71278895f09eaa3e4585ef06c40390b571c43403fdc3457f4c10618fe16389ac`), converted from RTF
to plain text with macOS `textutil -convert txt` and saved tab-separated with a renamed header;
no value was changed. The original header line is kept as the second comment line. The other
two supplementary files hold the individual core records (irregular ages, and resampled at
1 kyr) and are not used.

Columns (one row per kyr, ages 2–800 ka, 799 rows; HSO's "799 points"):

| Column | Original header | Meaning |
|---|---|---|
| `age_ka` | AGE (ka) | age, kyr before present |
| `n_cores` | NUMBER OF CORES | number of records contributing at that age |
| `sint800` | SINT-800 | the stacked relative palaeointensity |
| `stderr_2sigma` | STD. ERR. (2 SIGMA) | its standard error, 2σ |
| `stderr_1sigma` | STD. ERR. (1 SIGMA) | its standard error, 1σ |
| `vadm_1e22Am2` | VADM 10^22 Am^2 (Calib. 0-40 ka) | VADM, calibrated on the absolute data of the past 40 kyr, 10²² A m² |
| `vadm_stderr_1e22Am2` | STD ERR. VADM | its standard error (1σ), 10²² A m²; the s_i ± e_i of HSO p. 148 |

The record is unsigned (a relative palaeointensity stack carries no polarity); HSO's Fig. 1
flips the sign of the points older than the Brunhes–Matuyama reversal at 780 ka. Whole-record
check: mean VADM 5.875×10²² A m², dv = 0.0706 (HSO: 0.071 ± 0.003), median relative error 12%
(HSO: "as a rule e_i/s_i ∼ 0.1").

The data are the authors' and Nature's; they are reproduced here as published supplementary
material for the purpose of reproducing the analysis of HSO. Cite Guyodo & Valet (1999) when
using them.

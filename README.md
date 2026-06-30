# Brain Injury Mouse Model

This repository contains the computational models, analysis scripts, and experimental data used in the study:

**"Experimentally-Calibrated Poro-Hyperelastic Model of Brain Injury Reveals Astrocyte-Driven Edema Formation"**


## Overview

The repository provides all resources required to reproduce the computational simulations, model calibration, sensitivity analyses, and experimental data described in the manuscript.

## Repository Structure

```text
-- Abaqus-codes
-- Calibration-codes
-- Experimental-data
-- Plotting-codes
-- Sensitivity-codes
```

### Calibration-codes

Scripts used to calibrate model parameters against experimental measurements from mouse brain tissue sections (Figure 4 & 5): 

1. Copy the input file and the UEL from "Abaqus-codes" 
2. Change `CONFIGURATION_DATA["pop_size"]` and `CONFIGURATION_DATA["generations"]` to set the differential evolution run size.
3. Submit the bash file with: `qsub -pe smp <num of cpus> run_calibration.sh`
4. ROI handling.
   - `calibrate_de.py` expects `rois_random.json`.
   - If `rois_random.json` is missing, it is generated automatically before calibration.
5. Post-calibration sensitivity.
   - `sensitivity_de.py` evaluates sensitivity of the calibrated parameters.
   - It runs on the results stored in `calib_result_de.json`.
6. Calibration plotting.
   - Use `plot_calibration_results.py` to plot calibrated model outputs against the experimental data and objective function over runs.

### Experimental-data

Experimental datasets used for model calibration. These files support the imaging-based plots shown in Figures 6 and 7 of the paper.

### Plotting-codes

Scripts used to generate the figures and visualizations presented in the manuscript. These files produce Figure 3 (mesh study) and Figures 9–10 (parametric study) from the paper.

### Sensitivity-codes

Scripts for sensitivity analyses investigating the influence of injury location on the model by running Abaqus simulations at sampled center locations.

1. Copy the input file and the UEL from "Abaqus-codes" 
2. Submit the provided bash file with: `qsub run_sensitivity.sh`
3. Output is produced, including generated Abaqus ODB files and CSV exports to reproduce Figure 8 from the paper. 

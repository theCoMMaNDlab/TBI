# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import json
import os
import csv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

from calibrate_de import (
    CONFIGURATION_DATA,
    ensure_rois_file_exists,
    load_astrocytes_experimental_data,
    load_experimental_pressure_data,
    relative_squared_error,
    run_abaqus_and_extract_counts,
)

CUSTOM_PURD = [
    "#fef6f7", "#fde9ed", "#fcdbde", "#faccce", "#f7b9bf",
    "#f4a4b0", "#f08fa2", "#ec7894", "#e86187", "#e3497b",
    "#de3070", "#d71965", "#ce035b", "#b80353", "#a1034a",
    "#8a0341", "#720338", "#5b022e", "#440224", "#2d011a",
]

# ============================================================================
# Configuration
# ============================================================================
SENSITIVITY_CONFIG = {
    "calibrated_json": "calib_result_de.json",
    "scales": [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20],
    "out_dir": "sensitivity_outputs",
}

LABEL_MAP = {
    "k_s": r"$k_{\rm s}$",
    "Sa0": r"$S_{\rm {a0}}$",
    "delta_a": r"$\delta_{\rm a}$",
    "Sf0": r"$S_{\rm {f0}}$",
    "alpha_f": r"$\alpha_{\rm f}$",
    "epsilon": r"$\epsilon$",
    "beta": r"$\beta$",
}

# ============================================================================
# Core functions
# ============================================================================


def objective_for_params(params, param_names, input_template, uel_file, roi_file,
                         work_dir, times, exp_astro, exp_pressure, astrocyte_weight, pressure_weight):
    # Evaluate objective for a single parameter set (runs Abaqus and computes weighted error)
    astro_pred, pressure_pred = run_abaqus_and_extract_counts(params, param_names, input_template, uel_file, roi_file, work_dir, times)
    astro_error = relative_squared_error(exp_astro, astro_pred)
    pressure_error = relative_squared_error(exp_pressure, pressure_pred)
    return float(astrocyte_weight * astro_error + pressure_weight * pressure_error)


def build_single_parameter_sensitivity(best_params, param_names, scales, objective_eval):
    # Sweep each parameter over `scales` while holding others fixed and record objective
    results = []
    for pname in param_names:
        base = float(best_params[pname])
        tested_values = [float(base * s) for s in scales]
        losses = []
        for s in scales:
            params = dict(best_params)
            params[pname] = float(base * s)
            losses.append(float(objective_eval(params)))
        results.append({
            "parameter": pname,
            "tested_values": tested_values,
            "losses": losses,
        })
    return results


def save_single_parameter_csv(results, out_dir, scales):
    # Write sensitivity table to CSV
    out_path = os.path.join(out_dir, "single_parameter_sensitivity.csv")
    pct_labels = [f"{int(round(100 * s))}pct" for s in scales]
    headers = ["parameter"] + [f"value_{l}" for l in pct_labels] + [f"objective_{l}" for l in pct_labels]
    rows = [[item["parameter"]] + item["tested_values"] + item["losses"] for item in results]
    
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    return out_path


def plot_single_parameter_heatmap(results, out_path, scales):
    # Create and save heatmap of objective vs parameter perturbation
    mat = np.array([item["losses"] for item in results], dtype=float)
    params = [item["parameter"] for item in results]
    x_labels = [f"{int(round(100 * s))}%" for s in scales]
    # convert parameter names to LaTeX labels using LABEL_MAP
    y_labels = [LABEL_MAP.get(p, p) for p in params]

    fig, ax = plt.subplots(figsize=(4, max(2.0, 0.5 * len(params))))
    cmap = ListedColormap(CUSTOM_PURD)
    im = ax.imshow(mat, cmap=cmap, aspect="auto", interpolation="nearest")

    ax.set_xticks(np.arange(len(scales)))
    ax.set_xticklabels(x_labels, rotation=45, ha="center")
    ax.set_yticks(np.arange(len(params)))
    ax.set_yticklabels(y_labels)
    ax.tick_params(axis="both", which="major", labelsize=8)
    ax.set_xlabel("value change relative to best fit")
    ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=True, length=0)

    if mat.shape[0] > 1:
        ax.hlines(np.arange(0.5, mat.shape[0] - 0.5), -0.5, mat.shape[1] - 0.5,
                  colors="white", linewidth=1.5, zorder=3)

    cbar = fig.colorbar(im, ax=ax, shrink=0.9, fraction=0.12, pad=0.04, aspect=12)
    cbar.set_label(r"$f(\mathbf{\theta})$")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def main():
    out_dir = SENSITIVITY_CONFIG["out_dir"]
    scales = np.array(SENSITIVITY_CONFIG["scales"], dtype=float)
    os.makedirs(out_dir, exist_ok=True)

    with open(SENSITIVITY_CONFIG["calibrated_json"], "r") as f:
        payload = json.load(f)

    best_params = {k: float(v) for k, v in payload["best_params"].items()}
    param_names = CONFIGURATION_DATA["param_names"]
    times = np.array(CONFIGURATION_DATA["times"])

    ensure_rois_file_exists(CONFIGURATION_DATA["rois"])

    exp_astro = load_astrocytes_experimental_data(CONFIGURATION_DATA["exp_data"], times=times, group_name="TBI")
    exp_pressure, _ = load_experimental_pressure_data(CONFIGURATION_DATA["exp_pressure_data"], times=times)

    def objective_eval(params):
        return objective_for_params(
            params=params,
            param_names=param_names,
            input_template=CONFIGURATION_DATA["input_template"],
            uel_file=CONFIGURATION_DATA["uel_file"],
            roi_file=CONFIGURATION_DATA["rois"],
            work_dir=CONFIGURATION_DATA["work_dir"],
            times=times,
            exp_astro=exp_astro,
            exp_pressure=exp_pressure,
            astrocyte_weight=CONFIGURATION_DATA["astrocyte_weight"],
            pressure_weight=CONFIGURATION_DATA["pressure_weight"],
        )

    sensitivity_results = build_single_parameter_sensitivity(best_params, param_names, scales, objective_eval)
    csv_out = save_single_parameter_csv(sensitivity_results, out_dir, scales)
    fig_out = os.path.join(out_dir, "single_parameter_sensitivity_heatmap.png")
    plot_single_parameter_heatmap(sensitivity_results, fig_out, scales)

    summary = {
        "calibrated_json": SENSITIVITY_CONFIG["calibrated_json"],
        "scales": [float(s) for s in scales],
        "levels": int(len(scales)),
        "n_parameters": int(len(param_names)),
        "n_simulations": int(len(param_names) * len(scales)),
        "best_params": best_params,
        "objective_weights": {
            "astrocyte": CONFIGURATION_DATA["astrocyte_weight"],
            "pressure": CONFIGURATION_DATA["pressure_weight"],
        },
        "plot_file": fig_out,
        "csv_file": csv_out,
        "parameters": sensitivity_results,
    }

    with open(os.path.join(out_dir, "sensitivity_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

if __name__ == "__main__":
    main()

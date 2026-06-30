# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import json
import numpy as np
from scipy.optimize import differential_evolution
import time
import os
import csv
from visualize_rois import generate_random_rois
from abaqus_runs import (write_parameters_to_input, extract_average_data,
                         extract_per_roi_data, run_abaqus_job, save_per_roi_results)

# ============================================================================
# Configuration and Constants
# ============================================================================

# Optimization constraints
LARGE_ERROR_VALUE = 1e10  # Penalty value for failed simulations
EXP_ASTROCYTE = "experimental_astrocytes_density.csv"
EXP_PRESSURE = "experimental_pressure_mild.csv"
ASTROCYTE_WEIGHT = 3.0
PRESSURE_WEIGHT = 1.0
ROI_FILE = "rois_random.json"
PRESSURE_KEYWORD = "NT11"
ASTROCYTE_KEYWORD = "NT12"

# One-time ROI generation defaults (used only if ROI file does not exist)
INJURY_CENTER_X = 2000.0
INJURY_CENTER_Y = 4000.0
INNER_RADIUS = 0.0
OUTER_RADIUS = 3000.0
ANGLE_START = 0.0
ANGLE_END = 360.0
ROI_WIDTH = 1330.0
ROI_HEIGHT = 1330.0
N_ROIS = 3
DOMAIN_RADIUS = 5000.0

# Output files
CSV_OUTPUT = "nt12_per_roi_results.csv"
PRESSURE_CSV_OUTPUT = "nt11_per_roi_results.csv"
PARAM_HISTORY = "de_parameter_history.csv"
GEN_HISTORY = "de_generation_history.csv"
PROGRESS_LOG = "de_run_progress.log"

# ============================================================================
# Loading data
# ============================================================================

def ensure_rois_file_exists(roi_file_path):
    if os.path.exists(roi_file_path):
        return
    generate_random_rois(injury_center_x=INJURY_CENTER_X, injury_center_y=INJURY_CENTER_Y,
                         inner_radius=INNER_RADIUS, outer_radius=OUTER_RADIUS,
                         angle_start=ANGLE_START, angle_end=ANGLE_END, 
                         roi_width=ROI_WIDTH, roi_height=ROI_HEIGHT,
                         n_rois=N_ROIS, out_file=roi_file_path, radius=DOMAIN_RADIUS)

def load_astrocytes_experimental_data(csv_path, times, group_name):
    data = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        normalized = {name.strip().lower(): name for name in reader.fieldnames}
        group_col = normalized["group"]
        time_col = normalized["time_days"]
        value_col = normalized.get("astrocytes_density") or normalized["density"]

        for row in reader:
            grp = row[group_col].strip()
            time = float(row[time_col])
            value = float(row[value_col])
            data.setdefault((grp, time), []).append(value)

    means = {k: float(np.mean(v)) for k, v in data.items()}
    order_means = np.array([means[(group_name, float(t))] for t in times])
    return order_means


def load_experimental_pressure_data(csv_path, times):
    pressure_by_time = {}
    pressure_sd_by_time = {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            t_days = float(row["time_days"])
            pressure_by_time[t_days] = float(row["pressure"])
            pressure_sd_by_time[t_days] = float(row["SD"])

    exp_pressure = np.array([pressure_by_time[float(t)] for t in times], dtype=float)
    exp_pressure_sd = np.array([pressure_sd_by_time[float(t)] for t in times], dtype=float)
    return exp_pressure, exp_pressure_sd

def relative_squared_error(y_true, y_pred):
    errors = [(y_true[i] - y_pred[i]) / y_true[i] for i in range(len(y_true)) if y_true[i] != 0]
    if not errors:
        return LARGE_ERROR_VALUE
    return float(np.sum(np.array(errors) ** 2))

# ============================================================================
# Simulation execution
# ============================================================================

def run_abaqus_and_extract_counts(params, param_names, input_template, uel_file, roi_file, work_dir, times):
    if isinstance(params, (list, tuple, np.ndarray)):
        params = dict(zip(param_names, map(float, params)))

    input_out = os.path.splitext(input_template)[0] + "_de_run.inp"
    write_parameters_to_input(params, input_template, input_out)
    job_name = os.path.splitext(os.path.basename(input_out))[0]
    odb_path = run_abaqus_job(job_name, input_out, uel_file, work_dir=work_dir)
    with open(roi_file, "r") as f:
        rois = json.load(f)
    astrocyte_density = extract_average_data(odb_path, times, rois, ASTROCYTE_KEYWORD)
    pressure = extract_average_data(odb_path, times, rois, PRESSURE_KEYWORD)
    return astrocyte_density, pressure

# ============================================================================
# Optimization and result saving
# ============================================================================
def create_objective_function(param_names, input_template, uel_file, roi_file, work_dir, times,
                              exp_astrocyte_data, exp_pressure_data,
                              astrocyte_weight=ASTROCYTE_WEIGHT,
                              pressure_weight=PRESSURE_WEIGHT,
                              total_expected_runs=None,
                              progress_log_path=None):
    loss_history = []
    param_history = []
    run_counter = {"value": 0}

    def objective(x):
        run_counter["value"] += 1
        run_idx = run_counter["value"]
        if total_expected_runs is not None:
            run_label = f"run {run_idx}/{total_expected_runs}"
        else:
            run_label = f"run {run_idx}"

        param_history.append([float(v) for v in x])
        params_dict = dict(zip(param_names, x))
        astrocyte_pred, pressure_pred = run_abaqus_and_extract_counts(params_dict, param_names, input_template,
                                                                      uel_file, roi_file, work_dir, times)

        astrocyte_error = relative_squared_error(exp_astrocyte_data, astrocyte_pred)
        pressure_error = relative_squared_error(exp_pressure_data, pressure_pred)
        total_loss = astrocyte_weight * astrocyte_error + pressure_weight * pressure_error
        loss_history.append(float(total_loss))

        params_str = ", ".join(f"{name}={float(params_dict[name]):.6g}" for name in param_names)
        run_msg = f"{run_label} | loss={float(total_loss):.6g} | params: {params_str}"
        print(run_msg, flush=True)
        if progress_log_path:
            with open(progress_log_path, "a") as f:
                f.write(run_msg + "\n")

        return total_loss
    return objective, loss_history, param_history

def save_parameter_history(param_names, param_history, loss_history, output_file):
    n_rows = min(len(param_history), len(loss_history))
    headers = ["run", "loss"] + list(param_names)
    rows = [[idx, float(loss_history[idx])] + [float(v) for v in param_history[idx]]
            for idx in range(n_rows)]
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

def save_generation_history(param_names, generation_history, output_file):
    headers = ["generation", "evals_through_generation", "best_loss_so_far", "convergence"] + list(param_names)
    rows = [
        [
            int(row["generation"]),
            int(row["evals_through_generation"]),
            float(row["best_loss_so_far"]),
            float(row["convergence"]),
        ] + [float(v) for v in row["params"]]
        for row in generation_history
    ]
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

# ============================================================================
# Main per-ROI data extraction and saving functions
# ============================================================================

CONFIGURATION_DATA = {
    "param_names": ("k_s", "Sa0", "delta_a", "Sf0", "alpha_f", "epsilon", "beta"),
    "bounds": [(0.1, 10.0), (0.001, 0.1), (10, 100), (0.001, 0.1), (0.1, 10), (0.1, 5.0), (0.1, 5.0)],
    "times": [0.25, 1.0, 3.0, 7.0],
    "input_template": "2Dcircle.inp",
    "uel_file": "UEL_Poro_Cell.for",
    "work_dir": None,
    "pop_size": 7,
    "generations": 15,
    "seed": 42,
    "out": "calib_result_de.json",
    "exp_data": EXP_ASTROCYTE,
    "exp_pressure_data": EXP_PRESSURE,
    "astrocyte_weight": ASTROCYTE_WEIGHT,
    "pressure_weight": PRESSURE_WEIGHT,
    "rois": ROI_FILE,
}

def main():
    ensure_rois_file_exists(CONFIGURATION_DATA["rois"])

    start = time.time()
    param_names = CONFIGURATION_DATA["param_names"]
    bounds = np.array(CONFIGURATION_DATA["bounds"], dtype=float)
    times = np.array(CONFIGURATION_DATA["times"])

    exp_astrocyte_data = load_astrocytes_experimental_data(CONFIGURATION_DATA["exp_data"], times=times, group_name="TBI")
    exp_pressure_data, _ = load_experimental_pressure_data(CONFIGURATION_DATA["exp_pressure_data"], times=times)

    effective_n_params = int(np.sum(bounds[:, 0] != bounds[:, 1]))
    if effective_n_params <= 0:
        effective_n_params = len(param_names)
    evals_per_population = max(1, CONFIGURATION_DATA["pop_size"] * effective_n_params)
    total_runs_requested = max(1, CONFIGURATION_DATA["pop_size"] * CONFIGURATION_DATA["generations"] * effective_n_params)
    total_runs_scipy = max(1, (CONFIGURATION_DATA["generations"] + 1) * CONFIGURATION_DATA["pop_size"] * effective_n_params)

    print(f"Total runs (pop_size * generations * n_params): {total_runs_requested}", flush=True)
    print(f"Expected function evaluations in SciPy DE ((generations + 1) * pop_size * n_params): {total_runs_scipy}", flush=True)
    with open(PROGRESS_LOG, "w") as f:
        f.write(f"Total runs (pop_size * generations * n_params): {total_runs_requested}\n")
        f.write(
            "Expected function evaluations in SciPy DE "
            f"((generations + 1) * pop_size * n_params): {total_runs_scipy}\n"
        )

    objective_func, loss_history, param_history = create_objective_function(
        param_names, CONFIGURATION_DATA["input_template"], CONFIGURATION_DATA["uel_file"], CONFIGURATION_DATA["rois"],
        CONFIGURATION_DATA["work_dir"], times, exp_astrocyte_data, exp_pressure_data,
        astrocyte_weight=CONFIGURATION_DATA["astrocyte_weight"],
        pressure_weight=CONFIGURATION_DATA["pressure_weight"],
        total_expected_runs=total_runs_scipy,
        progress_log_path=PROGRESS_LOG,
    )

    generation_history = []
    generation_counter = {"value": 0}

    def de_callback(kk, convergence):
        generation_counter["value"] += 1
        generation_idx = generation_counter["value"]
        eval_cutoff = min(len(loss_history), (generation_idx + 1) * evals_per_population)
        best_loss_so_far = float(np.min(loss_history[:eval_cutoff])) if eval_cutoff > 0 else LARGE_ERROR_VALUE
        generation_history.append({
            "generation": generation_idx,
            "evals_through_generation": eval_cutoff,
            "best_loss_so_far": best_loss_so_far,
            "convergence": float(convergence),
            "params": [float(v) for v in kk],
        })
        return False

    result = differential_evolution(objective_func, bounds, popsize=CONFIGURATION_DATA["pop_size"],
                                    maxiter=CONFIGURATION_DATA["generations"], polish=False, seed=CONFIGURATION_DATA["seed"],
                                    callback=de_callback)

    save_parameter_history(param_names, param_history, loss_history, PARAM_HISTORY)
    save_generation_history(param_names, generation_history, GEN_HISTORY)
    best_run_index = int(np.argmin(loss_history)) if loss_history else None

    best_params = dict(zip(param_names, map(float, result.x)))
    astrocytes_sim_result, pressure_sim_result = run_abaqus_and_extract_counts(
        best_params, param_names, CONFIGURATION_DATA["input_template"], CONFIGURATION_DATA["uel_file"],
        CONFIGURATION_DATA["rois"], CONFIGURATION_DATA["work_dir"], times
    )

    astrocytes_sim_result = [float(x) for x in astrocytes_sim_result]
    pressure_sim_result = [float(x) for x in pressure_sim_result]
    input_out = os.path.splitext(CONFIGURATION_DATA["input_template"])[0] + "_de_run.inp"
    job_name = os.path.splitext(os.path.basename(input_out))[0]
    odb_path = os.path.join(CONFIGURATION_DATA["work_dir"] or os.getcwd(), f"{job_name}.odb")
    with open(CONFIGURATION_DATA["rois"], "r") as f:
        rois = json.load(f)

    per_roi_astrocyte_data = extract_per_roi_data(odb_path, None, rois, ASTROCYTE_KEYWORD)
    save_per_roi_results(per_roi_astrocyte_data, CSV_OUTPUT)

    per_roi_pressure_data = extract_per_roi_data(odb_path, None, rois, PRESSURE_KEYWORD)
    save_per_roi_results(per_roi_pressure_data, PRESSURE_CSV_OUTPUT)

    outcome = {
        "best_params": best_params,
        "function": float(result.fun),
        "num_iterations": int(result.nit),
        "message": result.message,
        "time_s": time.time() - start,
        "all_params": [float(x) for x in result.x],
        "astrocytes_sim_result": astrocytes_sim_result,
        "astrocytes_experimental_targets": [float(x) for x in exp_astrocyte_data],
        "pressure_sim_result": pressure_sim_result,
        "pressure_experimental_targets": [float(x) for x in exp_pressure_data],
        "objective_weights": {
            "astrocyte": float(CONFIGURATION_DATA["astrocyte_weight"]),
            "pressure": float(CONFIGURATION_DATA["pressure_weight"]),
        },
        "losses": [float(v) for v in loss_history],
        "best_run_index": best_run_index,
        "n_function_evals_logged": int(len(loss_history)),
        "n_generations_logged": int(len(generation_history)),
        "parameter_history": PARAM_HISTORY,
        "generation_history": GEN_HISTORY,
        "run_config": CONFIGURATION_DATA,
    }

    with open(CONFIGURATION_DATA["out"], "w") as f:
        json.dump(outcome, f, indent=2)

if __name__ == "__main__":
    main()

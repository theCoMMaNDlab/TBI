# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import os
import subprocess
import numpy as np
import time
import csv
from typing import Optional
from odbAccess import openOdb


# ============================================================================
# Input file parameter writing
# ============================================================================

def write_parameters_to_input(params, input_path, input_out_path):
    with open(input_path, "r") as f:
        lines = f.readlines()

    out_lines = []
    for line in lines:
        stripped = line.lstrip()
        for pname, val in params.items():
            if stripped.startswith(pname):
                left = line.split("=", 1)[0]
                out_lines.append(f"{left.rstrip()} = {val}\n")
                break
        else:
            out_lines.append(line)

    with open(input_out_path, "w") as f:
        f.writelines(out_lines)


# ============================================================================
# ABAQUS JOB EXECUTION
# ============================================================================

def run_abaqus_job(job_name, input_file, uel_file, work_dir=None, timeout=None,
                   poll_interval=5.0, cpus=int(os.environ.get("NSLOTS", 1))):
    
    cwd = work_dir or os.getcwd()
    command = (
        f"abaqus double job={job_name} input={input_file} user={uel_file}"
        f" cpus={cpus} mp_mode=threads"
    )
    remove_exts = [".odb", ".dat", ".msg", ".sta", ".prt", ".com", ".lck", ".log", ".res"]
    for ext in remove_exts:
        path = os.path.join(cwd, job_name + ext)
        if os.path.exists(path):
            os.remove(path)
    subprocess.run(command, shell=True, check=True, cwd=cwd, timeout=timeout)
    odb_path = os.path.join(cwd, f"{job_name}.odb")
    odb_dir = os.path.dirname(odb_path)
    odb_base = os.path.splitext(os.path.basename(odb_path))[0]
    lck_path = os.path.join(odb_dir, odb_base + ".lck")
    while True:
        if os.path.exists(lck_path):
            time.sleep(poll_interval)
            continue
        if os.path.exists(odb_path):
            return odb_path
        time.sleep(poll_interval)


# ============================================================================
# ODB results extraction
# ============================================================================

def extract_average_data(odb_path, times, rois, output_name):
    odb = openOdb(path=odb_path)
    root = odb.rootAssembly
    instance = root.instances[list(root.instances.keys())[0]]
    node_coords = {
        node.label: (node.coordinates[0], node.coordinates[1])
        for node in instance.nodes
    }

    all_frames = [frame for step in odb.steps.values() for frame in step.frames]

    preds = []
    for target_time in times:
        closest_frame = min(all_frames, key=lambda f: abs(float(f.frameValue) - target_time))
        node_values = {
            value.nodeLabel: float(value.data)
            for value in closest_frame.fieldOutputs[output_name].values
            if getattr(value, "nodeLabel", None) is not None
        }

        roi_means = []
        for (cx, cy, width_m, height_m) in rois:
            x_min, x_max = cx - width_m / 2.0, cx + width_m / 2.0
            y_min, y_max = cy - height_m / 2.0, cy + height_m / 2.0
            vals = [node_values[lbl] for lbl, (x, y) in node_coords.items()
                    if x_min <= x <= x_max and y_min <= y <= y_max and lbl in node_values]
            roi_means.append(float(np.mean(vals)) if vals else 0.0)

        preds.append(float(np.mean(roi_means)) if roi_means else 0.0)

    odb.close()
    return np.array(preds)

def extract_per_roi_data(odb_path, times, rois, output_name):
    odb = openOdb(path=odb_path)
    instance = odb.rootAssembly.instances[list(odb.rootAssembly.instances.keys())[0]]
    node_coords = {node.label: (node.coordinates[0], node.coordinates[1])
                   for node in instance.nodes}

    all_frames = [frame for step in odb.steps.values() for frame in step.frames]

    def get_roi_values(frame):
        node_values = {
            value.nodeLabel: float(value.data)
            for value in frame.fieldOutputs[output_name].values
            if getattr(value, "nodeLabel", None) is not None
        }
        roi_means = []
        for (cx, cy, width_m, height_m) in rois:
            x_min, x_max = cx - width_m / 2.0, cx + width_m / 2.0
            y_min, y_max = cy - height_m / 2.0, cy + height_m / 2.0
            vals = [node_values[lbl] for lbl, (x, y) in node_coords.items()
                    if x_min <= x <= x_max and y_min <= y <= y_max and lbl in node_values]
            roi_means.append(float(np.mean(vals)) if vals else 0.0)
        return roi_means

    per_roi_data = {}
    if times is None:
        for frame in all_frames:
            per_roi_data[float(frame.frameValue)] = get_roi_values(frame)
    else:
        for target_time in times:
            closest_frame = min(all_frames, key=lambda f: abs(float(f.frameValue) - target_time))
            per_roi_data[float(target_time)] = get_roi_values(closest_frame)

    odb.close()
    return per_roi_data

def save_per_roi_results(per_roi_data, output_file):
    sorted_times = sorted(per_roi_data.keys())
    num_rois = len(per_roi_data[sorted_times[0]])
    headers = ["time"] + [f"ROI{i+1}_avg" for i in range(num_rois)] + ["average_all_rois"]
    rows = [
        [t] + per_roi_data[t] + [float(np.mean(per_roi_data[t]))]
        for t in sorted_times
    ]
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)

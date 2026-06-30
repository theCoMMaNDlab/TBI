# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import csv
import math
import os
import random
import re
import subprocess
import time
from odbAccess import openOdb

#---------------------------------------------------------
# Configuration
#---------------------------------------------------------
OUTPUT_POINTS = (150.0, 300.0, 450.0)
MAX_DISTANCE = 300.0
N_SAMPLES = 20
INSTANCE_NAME = 'PART-1-1'
STEP_NAME = 'Swell'
UEL_FILE = 'UEL_Poro_Cell.for'
INPUT_FILE = '2Dcircle.inp'
RANDOM_SEED = 42
#---------------------------------------------------------
# Abaqus functions for sensitivity analysis
#---------------------------------------------------------
def write_new_center_inputs(names, centers, input_path, out_dir, start_index):
    input_path = str(input_path)
    lines = open(input_path, 'r').readlines()

    params = {}
    if names:
        names_escaped = [re.escape(n) for n in names]
        pattern = re.compile(
            r'^\s*(?P<name>' + '|'.join(names_escaped) + r')\b\s*=\s*(?P<value>[^#\n\r]*)'
        )
        for line in lines:
            stripped = line.lstrip()
            if not stripped or stripped[0] in ('!', '*') or stripped.startswith('//'):
                continue
            m = pattern.match(line)
            if not m:
                continue
            name = m.group('name')
            raw_val = m.group('value').split('!')[0].split('//')[0].strip()
            val = float(raw_val) if raw_val.replace('.', '', 1).lstrip('-').isdigit() else raw_val.strip().strip('"\'')
            params[name] = val

    out_paths = []
    if centers:
        out_dir = str(out_dir)
        os.makedirs(out_dir, exist_ok=True)
        input_stem = os.path.splitext(os.path.basename(input_path))[0]
        for i, (xc, yc) in enumerate(centers, start_index):
            center_params = {'x_center': xc, 'y_center': yc}
            out_lines = []
            for line in lines:
                replaced_line = line
                for pname, pval in center_params.items():
                    m = re.match(
                        r'^(?P<indent>\s*)(?P<name>' + re.escape(pname) + r')\b\s*=\s*(?P<val>[^#\n\r]*)(?P<rest>.*)$',
                        line,
                    )
                    if m:
                        indent = m.group('indent') or ''
                        rest = m.group('rest') or ''
                        replaced_line = f"{indent}{pname} = {pval}{rest}\n"
                        break
                out_lines.append(replaced_line)
            out_path = os.path.join(out_dir, f"{input_stem}_center_{i}.inp")
            open(out_path, 'w').writelines(out_lines)
            out_paths.append(out_path)

    return params, out_paths


def run_abaqus_job(job_name, input_file, uel_file, work_dir=None, timeout=None,
                   poll_interval=5.0, cpus=int(os.environ.get("NSLOTS", 1))):
    cwd = work_dir or os.getcwd()

    command = (
        f"abaqus double job={job_name} input={input_file} user={uel_file}"
        f" cpus={cpus} mp_mode=threads"
    )

    for ext in [".odb", ".dat", ".msg", ".sta", ".prt", ".com", ".lck", ".log", ".res"]:
        path = os.path.join(cwd, job_name + ext)
        if os.path.exists(path):
            os.remove(path)

    subprocess.run(command, shell=True, check=True, cwd=cwd, timeout=timeout)

    odb_path = os.path.join(cwd, f"{job_name}.odb")
    lck_path = os.path.join(cwd, f"{job_name}.lck")

    while True:
        if os.path.exists(lck_path):
            time.sleep(poll_interval)
            continue
        if os.path.exists(odb_path):
            return odb_path
        time.sleep(poll_interval)

def extract_outputs_from_odb(odb_path, instance_name, step_name, baseline_center,
                             sample_radial_offsets=OUTPUT_POINTS):
    csv_path = f"{odb_path}.csv"
    odb = openOdb(path=odb_path)
    instance = odb.rootAssembly.instances[instance_name]
    bx, by = baseline_center
    angle = math.pi * 5 / 4

    sample_points = [(bx, by)] + [
        (bx + r * math.cos(angle), by + r * math.sin(angle))
        for r in sample_radial_offsets
    ]

    node_positions = [
        (node.label, node.coordinates[0], node.coordinates[1] if len(node.coordinates) > 1 else 0.0)
        for node in instance.nodes
    ]

    closest_nodes = [
        (min(node_positions, key=lambda nd: math.hypot(nd[1] - tx, nd[2] - ty))[0], tx, ty)
        for tx, ty in sample_points
    ]

    rows = []
    for frame in odb.steps[step_name].frames:
        nt11_map = {v.nodeLabel: v.data for v in frame.fieldOutputs['NT11'].values}
        nt12_map = {v.nodeLabel: v.data for v in frame.fieldOutputs['NT12'].values}
        u_map = {v.nodeLabel: math.sqrt(v.data[0] ** 2 + v.data[1] ** 2) for v in frame.fieldOutputs['U'].values}
        row = [frame.frameValue]
        for node_label, _, _ in closest_nodes:
            row.extend([nt11_map.get(node_label), nt12_map.get(node_label), u_map.get(node_label)])
        rows.append(row)

    headers = ['time']
    for _, px, py in closest_nodes:
        point_label = f"({px:.1f},{py:.1f})"
        headers.extend([f"{point_label}_NT11", f"{point_label}_NT12", f"{point_label}_UMAG"])

    writer = csv.writer(open(csv_path, 'w', newline=''))
    writer.writerow(headers)
    writer.writerows(rows)

    odb.close()
    return csv_path

def run_center_simulation(center, sample_index, input_file, uel_path, baseline_center):
    input_file = str(input_file)
    input_dir = os.path.dirname(input_file) or '.'
    _, input_created = write_new_center_inputs([], [center], input_file, input_dir, sample_index)
    inp_path = input_created[0]
    job_name = os.path.splitext(os.path.basename(inp_path))[0]
    odb_path = run_abaqus_job(job_name, inp_path, uel_path, work_dir=input_dir)
    extract_outputs_from_odb(odb_path, INSTANCE_NAME, STEP_NAME, baseline_center=baseline_center)

#---------------------------------------------------------
# Main
#---------------------------------------------------------

def main():
    parameters, _ = write_new_center_inputs(['x_center', 'y_center'], None, INPUT_FILE, None, 0)
    baseline_center = (float(parameters['x_center']), float(parameters['y_center']))
    bx, by = baseline_center

    random.seed(RANDOM_SEED)
    pts = []
    for _ in range(N_SAMPLES):
        radius = MAX_DISTANCE * math.sqrt(random.random())
        theta = random.uniform(0, 2 * math.pi)
        pts.append((bx + radius * math.cos(theta), by + radius * math.sin(theta)))

    input_dir = os.path.dirname(INPUT_FILE) or '.'
    writer = csv.writer(open(os.path.join(input_dir, 'sample_centers.csv'), 'w', newline=''))
    writer.writerow(['sample', 'x_coord', 'y_coord'])
    writer.writerow(['baseline', bx, by])
    writer.writerows([[idx, px, py] for idx, (px, py) in enumerate(pts, 1)])

    run_center_simulation(baseline_center, sample_index=0, input_file=INPUT_FILE,
                          uel_path=UEL_FILE, baseline_center=baseline_center)
    for idx, center in enumerate(pts, 1):
        run_center_simulation(center, sample_index=idx, input_file=INPUT_FILE,
                              uel_path=UEL_FILE, baseline_center=baseline_center)


# ============================================================================
if __name__ == '__main__':
    main()

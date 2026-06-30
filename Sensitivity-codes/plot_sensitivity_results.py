# Author: Fatemeh Jafarabadi, Apr 2026
# ============================================================================
import csv
import glob
import math
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
#---------------------------------------------------------
# Configuration
#---------------------------------------------------------
BASELINE_CENTER = (2000.0, 4000.0)
MAX_DISTANCE = 300.0

#---------------------------------------------------------
# Functions
#---------------------------------------------------------
def read_data_multipoints(csv_output):
    reader = csv.reader(open(csv_output, newline=''))
    header = next(reader)

    locations = [
        (float(m.group(1)), float(m.group(2)))
        for col in header[1:]
        if (m := re.match(r'\(([^,]+),([^)]+)\)', col.strip('"')))
    ]

    t, data_by_location = [], [[] for _ in range(len(locations))]
    for row in reader:
        t.append(float(row[0]))
        for i, loc_data in enumerate(data_by_location):
            loc_data.append(float(row[i + 1]))

    return np.array(t), locations, [np.array(d) for d in data_by_location]

def read_combined_outputs(csv_output):
    reader = csv.reader(open(csv_output, newline=''))
    header = next(reader)

    col_info = []
    locs_by_field = {'nt11': [], 'nt12': [], 'umag': []}
    for col in header[1:]:
        m = re.match(r'^\(([^,]+),([^)]+)\)_(NT11|NT12|UMAG)$', col.strip('"'))
        if not m:
            continue
        key = m.group(3).lower()
        locs_by_field[key].append((float(m.group(1)), float(m.group(2))))
        col_info.append((key, len(locs_by_field[key]) - 1))

    data_by_field = {k: [[] for _ in v] for k, v in locs_by_field.items()}
    t = []
    for row in reader:
        t.append(float(row[0]))
        for col_idx, (key, loc_idx) in enumerate(col_info, start=1):
            data_by_field[key][loc_idx].append(float(row[col_idx]))

    t = np.array(t)
    return {
        't_nt11': t, 'locations_nt11': locs_by_field['nt11'], 'data_nt11': [np.array(d) for d in data_by_field['nt11']],
        't_nt12': t, 'locations_nt12': locs_by_field['nt12'], 'data_nt12': [np.array(d) for d in data_by_field['nt12']],
        't_umag': t, 'locations_umag': locs_by_field['umag'], 'data_umag': [np.array(d) for d in data_by_field['umag']],
    }

def compute_mape_by_location(baseline_run, all_runs, field_type):
    t_key, data_key = f't_{field_type}', f'data_{field_type}'
    base_t, base_data = baseline_run[t_key], baseline_run[data_key]
    eps = 1e-12

    mape_by_location = [[] for _ in range(len(base_data))]
    for run in all_runs:
        for loc_idx, loc_data in enumerate(base_data):
            run_interp = np.interp(base_t, run[t_key], run[data_key][loc_idx])
            denominator = np.where(np.abs(loc_data) < eps, eps, loc_data)
            mape_by_location[loc_idx].append(np.mean(np.abs((run_interp - loc_data) / denominator) * 100.0))

    return mape_by_location


def create_boxplot(locations, mape_by_location, colors, field_type):
    fig, ax = plt.subplots(figsize=(6, 6))
    loc_distances = [math.sqrt((x - BASELINE_CENTER[0])**2 + (y - BASELINE_CENTER[1])**2) for x, y in locations]
    labels = [f'{d / MAX_DISTANCE:.1f}$r_{{\\rm s}}$' for d in loc_distances]

    ax.boxplot(mape_by_location, tick_labels=labels, patch_artist=True, widths=0.6,
               boxprops=dict(facecolor='lightgray', color='black'),
               medianprops=dict(color='black', linewidth=2),
               whiskerprops=dict(color='black'),
               capprops=dict(color='black'),
               whis=[0, 100], showfliers=False)

    rng = np.random.RandomState(1)
    for loc_idx in range(len(locations)):
        jitter = (rng.rand(len(mape_by_location[loc_idx])) - 0.5) * 0.18
        ax.scatter((loc_idx + 1) + jitter, mape_by_location[loc_idx],
                   c=colors, s=75, edgecolors='k', linewidths=0.5, zorder=10)

    sm = plt.cm.ScalarMappable(cmap=plt.cm.Blues, norm=plt.Normalize(0, MAX_DISTANCE))
    cbar = fig.colorbar(sm, ax=ax, orientation='horizontal', pad=0.12, fraction=0.08)
    cbar.set_ticks([0, 100, 200, 300, 400, 500])
    cbar.set_label('Distance from Baseline', fontsize=12)
    cbar.ax.tick_params(labelsize=12)

    ax.set_ylabel('MAPE (%)', fontsize=14)
    ax.set_ylim(-0.05, 3.05)
    ax.set_yticks([0.0, 0.50, 1.0, 1.5, 2.0, 2.5, 3.0])
    ax.tick_params(axis='both', labelsize=14)
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis='y', which='minor', length=3)
    ax.legend(handles=[plt.Line2D([], [], marker='None', linestyle='None')],
              loc='upper left', frameon=False, fontsize=16, prop={'weight': 'bold'})

    plt.tight_layout()
    plt.savefig(f'sensitivity_{field_type.upper()}.png', dpi=150)
    plt.show()

#---------------------------------------------------------
# Main
#---------------------------------------------------------
def main():
    fields_outputs = sorted(glob.glob('*_center_*.odb.csv'))
    sample_centers = {}
    with open('sample_centers.csv', newline='') as f:
        for row in csv.DictReader(f):
            sample_key = row['sample']
            center = (float(row['x_coord']), float(row['y_coord']))
            sample_centers[sample_key] = center
            if sample_key == 'baseline':
                sample_centers['0'] = center

    runs_data = []
    for output_files in fields_outputs:
        sim_num = re.search(r'center_(\d+)', output_files)
        sample_id = sim_num.group(1) if sim_num else None
        cx, cy = sample_centers.get(sample_id, (math.nan, math.nan))
        dist = math.sqrt((cx - BASELINE_CENTER[0])**2 + (cy - BASELINE_CENTER[1])**2) if not math.isnan(cx) else math.nan

        fields = read_combined_outputs(output_files)
        runs_data.append({'file_fields': output_files, 'center': (cx, cy), 'distance': dist, **fields})

    baseline_run = next((r for r in runs_data if re.search(r'_center_0\.odb\.csv$', r['file_fields'])), None)
    analysis_runs = [r for r in runs_data if r is not baseline_run]
    distances = np.array([r['distance'] if not math.isnan(r['distance']) else 0 for r in analysis_runs])
    colors = plt.cm.Blues(plt.Normalize(0, MAX_DISTANCE)(distances))

    for field_type in ['nt11', 'nt12', 'umag']:
        mape_by_location = compute_mape_by_location(baseline_run, analysis_runs, field_type)
        create_boxplot(baseline_run[f'locations_{field_type}'], mape_by_location, colors, field_type)

if __name__ == '__main__':
    main()

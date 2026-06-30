# Author: Fatemeh Jafarabadi, Apr 2026
# ============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy import stats as sp_stats

# ------------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------------

if_df = pd.read_csv('IF_data.csv')
if_df['group'] = if_df['group'].apply(lambda c: 'TBI' if 'TBI' in str(c) else 'Control')
if_df['cell_density'] = if_df['density']

mri_df = pd.read_csv('MRI_data.csv').dropna(subset=['time_days', 'group'])
mri_df['group'] = mri_df['group'].apply(lambda c: 'TBI' if 'TBI' in str(c) else 'Control')
mri_df = mri_df[mri_df['Segment'] == 'Ratio'].copy()
mri_df['Surface mm2 (CS)'] = pd.to_numeric(mri_df['Surface mm2 (CS)'], errors='coerce')
mri_df = mri_df[mri_df['time_days'].isin([0.25, 1, 3, 7])]
mri_df['time_days'] = mri_df['time_days']

# ------------------------------------------------------------------
# Functions 
# ------------------------------------------------------------------

def p_to_stars(p_value):
    if not np.isfinite(p_value):
        return 'n.s.'
    if p_value < 0.001:
        return '***'
    if p_value < 0.01:
        return '**'
    if p_value < 0.05:
        return '*'
    return 'n.s.'


def plot_grouped_comparison(df, *, value_col, output_path, ylabel, scale=1.0, y_ticks=None, scale_text=None):
    group_col = 'group'
    time_col = 'time_days'
    groups = ['Control', 'TBI']
    labels = ['Control', 'TBI']
    colors = ['deepskyblue', "#b80354d9"]
    time_points = [0.25, 1, 3, 7]
    
    stats = df.groupby([group_col, time_col])[value_col].agg(['mean', 'sem']).reset_index()
    times = list(time_points) if time_points is not None else sorted(stats[time_col].dropna().unique())
    x = np.arange(len(times))
    group_width = 0.5
    bar_width = group_width / max(1, len(groups))
    offsets = np.linspace(-group_width / 2 + bar_width / 2, group_width / 2 - bar_width / 2, len(groups)) if len(groups) > 1 else np.array([0.0])

    series = {}
    for group in groups:
        group_stats = stats[stats[group_col] == group]
        means = []
        sems = []
        for time_point in times:
            row = group_stats[group_stats[time_col] == time_point]
            means.append(row['mean'].values[0] * scale if len(row) else 0.0)
            sems.append(row['sem'].values[0] * scale if len(row) else 0.0)
        series[group] = {'means': means, 'sems': sems}

    fig, axis = plt.subplots(figsize=(8, 6))
    for idx, group in enumerate(groups):
        axis.bar(x + offsets[idx], series[group]['means'], bar_width, yerr=series[group]['sems'],
                 label=labels[idx], color=colors[idx], alpha=0.5, edgecolor=colors[idx],
                 linewidth=2.0, capsize=7.0, error_kw={'elinewidth': 1.5, 'capthick': 1.5}, zorder=3)

    for i, time_point in enumerate(times):
        group_a = df[(df[group_col] == groups[0]) & (df[time_col] == time_point)][value_col].dropna().values
        group_b = df[(df[group_col] == groups[1]) & (df[time_col] == time_point)][value_col].dropna().values
        if len(group_a) > 0 and len(group_b) > 0:
            _, p_value = sp_stats.ttest_ind(group_a, group_b)
            label = p_to_stars(p_value)
            top_a = series[groups[0]]['means'][i] + series[groups[0]]['sems'][i]
            top_b = series[groups[1]]['means'][i] + series[groups[1]]['sems'][i]
            y_pos = max(top_a, top_b) + 0.3
            left = x[i] + offsets[0]
            right = x[i] + offsets[1]
            axis.plot([left, left, right, right], [y_pos, y_pos + 0.5, y_pos + 0.5, y_pos], 'k-', linewidth=2.0, zorder=5)
            axis.text(x[i], y_pos + 0.5 + 0.03, label, ha='center', va='bottom',
                          fontsize=12 if label == 'n.s.' else 16, zorder=6, fontweight='bold')

    axis.set_xlabel("days (post-injury)", fontsize=16)
    axis.set_ylabel(ylabel, fontsize=16)
    axis.set_xticks(x)
    axis.set_xticklabels([str(time_point) for time_point in times], fontsize=16)
    axis.set_yticks(y_ticks)
    axis.set_yticklabels([str(tick) for tick in y_ticks], fontsize=16)
    axis.yaxis.set_minor_locator(AutoMinorLocator(5))
    axis.tick_params(axis='y', which='minor', length=2)
    axis.text(-0.01, 1.02, scale_text, transform=axis.transAxes, fontsize=14)
    max_height = max((max(np.array(series[group]['means']) + np.array(series[group]['sems'])) if series[group]['means'] else 0.0) for group in groups)
    axis.set_ylim(0, max_height * 1.5 if max_height > 0 else 1)
    axis.legend(fontsize=14, loc='upper right')
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)

    return axis


plot_grouped_comparison(if_df, value_col='cell_density', output_path='IF.png', ylabel=r'$\bar{a}_{\rm exp}\  [1/{\mu\mathrm{m}}^2]$',
                        scale=1e4, y_ticks=[0.0, 2.0, 4.0, 6.0, 8.0], scale_text=r'(10$^{-4}$)')


plot_grouped_comparison(mri_df, value_col='Surface mm2 (CS)', output_path='MRI.png', ylabel='VR %',
                        scale=100, y_ticks=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0])


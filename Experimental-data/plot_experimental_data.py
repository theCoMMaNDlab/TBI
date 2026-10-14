
# Author: Fatemeh Jafarabadi, Apr 2026
# Updated: Fatemeh Jafarabadi, Oct 2026
# ============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from scipy import stats as sp_stats
from statsmodels.stats.multitest import multipletests


# ------------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------------

if_df = pd.read_csv('IF_data.csv')
if_df['group'] = if_df['group'].apply(lambda c: 'TBI' if 'TBI' in str(c).upper() else 'Control')
if_df['cell_density'] = pd.to_numeric(if_df['density'], errors='coerce')
if_df['animal_id'] = if_df['experiment'].astype(str) + '_' + if_df['ID'].str.split('-').str[0]

mri_df = pd.read_csv('MRI_data.csv').dropna(subset=['time_days', 'group'])
mri_df['group'] = mri_df['group'].apply(lambda c: 'TBI' if 'TBI' in str(c).upper() else 'Control')
mri_df = mri_df[mri_df['Segment'] == 'Ratio'].copy()
mri_df['Surface mm2 (CS)'] = pd.to_numeric(mri_df['Surface mm2 (CS)'], errors='coerce')
mri_df['time_days'] = pd.to_numeric(mri_df['time_days'], errors='coerce')
mri_df = mri_df[mri_df['time_days'].isin([0.25, 1, 3, 7])].dropna(subset=['Surface mm2 (CS)'])

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


def nested_total_sd(df, group_col, time_col, value_col):
    out = []
    for (g, t), sub in df.groupby([group_col, time_col]):
        per = sub.groupby('animal_id')[value_col].agg(['mean', 'var', 'count'])
        per = per[per['count'] > 0]
        n = per['count'].to_numpy(float)
        ok = n >= 2
        if ok.sum() == 0 or len(per) < 2:
            out.append((g, t, np.nan)); continue

        # pooled within-animal variance (weights n_i - 1)
        s2_w = np.sum((n[ok] - 1) * per['var'].to_numpy()[ok]) / np.sum(n[ok] - 1)

        # between-animal variance, corrected for within-animal noise in the means
        var_means = per['mean'].var(ddof=1)
        s2_b = max(var_means - np.mean(s2_w / n), 0.0)

        out.append((g, t, np.sqrt(s2_b + s2_w)))
    return pd.DataFrame(out, columns=[group_col, time_col, 'sd_nested'])


def plot_grouped_comparison(df, *, value_col, output_path, ylabel, scale=1.0, y_ticks=None, scale_text=None):
    group_col = 'group'
    time_col = 'time_days'
    groups = ['Control', 'TBI']
    labels = ['Control', 'TBI']
    colors = ['deepskyblue', '#b80354d9']
    times = [0.25, 1, 3, 7]

    # --------------------------------------------------------------
    # Define the experimental unit for descriptive statistics.
    # IF: average ROI values within each animal first.
    # MRI: each ratio row represents one animal at one time point.
    # --------------------------------------------------------------
    if 'animal_id' in df.columns:
        animal_df = df.groupby([group_col, time_col, 'animal_id'], as_index=False)[value_col].mean()
    else:
        animal_df = df.copy()

    # Calculate group mean, SD, sample size across animals
    stats_df = animal_df.groupby([group_col, time_col])[value_col].agg(['mean', 'std', 'count']).reset_index()

    if 'animal_id' in df.columns:
        nested = nested_total_sd(df, group_col, time_col, value_col)
        stats_df = stats_df.merge(nested, on=[group_col, time_col], how='left')
        stats_df['sd'] = stats_df['sd_nested']
    else:
        stats_df['sd'] = stats_df['std']   # MRI: one value per animal, nothing to propagate
    x = np.arange(len(times))
    group_width = 0.5
    bar_width = group_width / len(groups)
    offsets = np.linspace(
        -group_width / 2 + bar_width / 2,
        group_width / 2 - bar_width / 2,
        len(groups)
    )

    series = {}
    for group in groups:
        group_stats = stats_df[stats_df[group_col] == group]
        means, sd = [], []

        for time_point in times:
            row = group_stats[group_stats[time_col] == time_point]

            if len(row):
                means.append(row['mean'].iloc[0] * scale)
                sd_value = row['sd'].iloc[0]
                sd.append(
                    sd_value * scale if np.isfinite(sd_value) else 0.0
                )
            else:
                means.append(0.0)
                sd.append(0.0)

        series[group] = {'means': means, 'sd': sd}

    # --------------------------------------------------------------
    # Statistical tests: Welch's independent two-sample t-test
    # at each time point, followed by Holm correction across time.
    # --------------------------------------------------------------
    raw_p_values = []

    for time_point in times:
        subset = animal_df[animal_df[time_col] == time_point]
        group_a = subset.loc[subset[group_col] == groups[0], value_col].dropna().to_numpy()
        group_b = subset.loc[subset[group_col] == groups[1], value_col].dropna().to_numpy()

        if len(group_a) >= 2 and len(group_b) >= 2:
            _, p_value = sp_stats.ttest_ind(
                group_a, group_b, equal_var=False
            )
        else:
            p_value = np.nan

        raw_p_values.append(p_value)

    adjusted_p_values = np.full(len(raw_p_values), np.nan)
    valid = np.isfinite(raw_p_values)

    if valid.any():
        adjusted_p_values[valid] = multipletests(
            np.asarray(raw_p_values)[valid],
            method='holm'
        )[1]

    # --------------------------------------------------------------
    # Plot bars and SD error bars
    # --------------------------------------------------------------
    fig, axis = plt.subplots(figsize=(8, 6))

    for idx, group in enumerate(groups):
        axis.bar(x + offsets[idx], series[group]['means'], bar_width, yerr=series[group]['sd'], label=labels[idx],
            color=colors[idx], alpha=0.5, edgecolor=colors[idx], linewidth=2.0, capsize=7.0, error_kw={'elinewidth': 1.5, 'capthick': 1.5}, zorder=3)

    # Add significance brackets using Holm-adjusted p-values
    for i, time_point in enumerate(times):
        label = p_to_stars(adjusted_p_values[i])

        top_a = (series[groups[0]]['means'][i] + series[groups[0]]['sd'][i])
        top_b = (series[groups[1]]['means'][i] + series[groups[1]]['sd'][i])

        y_pos = max(top_a, top_b) + 0.3
        left = x[i] + offsets[0]
        right = x[i] + offsets[1]

        axis.plot([left, left, right, right], [y_pos, y_pos + 0.5, y_pos + 0.5, y_pos], 'k-', linewidth=2.0, zorder=5)
        axis.text(x[i], y_pos + 0.53, label, ha='center', va='bottom',
                  fontsize=12 if label == 'n.s.' else 16, zorder=6, fontweight='bold')

    axis.set_xlabel("days (post-injury)", fontsize=18)
    axis.set_ylabel(ylabel, fontsize=18)
    axis.set_xticks(x)
    axis.set_xticklabels([str(t) for t in times], fontsize=16)

    axis.set_yticks(y_ticks)
    axis.set_yticklabels([str(tick) for tick in y_ticks], fontsize=16)

    axis.yaxis.set_minor_locator(AutoMinorLocator(5))
    axis.tick_params(axis='y', which='minor', length=2)

    if scale_text is not None:
        axis.text(-0.01, 1.02, scale_text, transform=axis.transAxes, fontsize=14)

    axis.legend(fontsize=16, loc='upper right')
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

    results = pd.DataFrame({
        'time_days': times,
        'raw_p': raw_p_values,
        'holm_adjusted_p': adjusted_p_values,
        'significance': [p_to_stars(p) for p in adjusted_p_values]
    })
    print(f'\nStatistical results: {output_path}')
    print(results.to_string(index=False))

    return results


# ------------------------------------------------------------------
# Generate figures
# ------------------------------------------------------------------

if_results = plot_grouped_comparison(if_df, value_col='cell_density', output_path='IF.png', ylabel=r'$\bar{a}_{\rm exp}\  [1/{\mu\mathrm{m}}^2]$',
                                      scale=1e4, y_ticks=[0.0, 2.0, 4.0, 6.0, 8.0], scale_text=r'(10$^{-4}$)')

mri_results = plot_grouped_comparison( mri_df, value_col='Surface mm2 (CS)', output_path='MRI.png', ylabel='VR (%)',
                                       scale=100, y_ticks=[0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
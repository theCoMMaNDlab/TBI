# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import json
from pathlib import Path

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================================
# CONFIGURATION
# ============================================================================
DENSITY_SCALE = 10000.0  
SIM_COLOR = "#164f7f"  
EXP_COLOR = "firebrick"  

def plot_density(base):
	sim_file = os.path.join(base, "nt12_per_roi_results.csv")
	exp_file = os.path.join(base, "experimental_astrocytes_density.csv")

	# ========================================================================
	# Load Simulation Data
	# ========================================================================
	df_sim = pd.read_csv(sim_file)
	time_sim = df_sim["time"].values

	roi_columns = [col for col in df_sim.columns if col.startswith("ROI") and col.endswith("_avg")]
	roi_data = df_sim[roi_columns].values

	mean_sim = np.mean(roi_data, axis=1)
	std_sim = np.std(roi_data, axis=1)

	# ========================================================================
	# Load and Process Experimental Data
	# ========================================================================
	df_exp = pd.read_csv(exp_file)
	df_exp = df_exp[df_exp["group"] == "TBI"]
	df_exp_stats = df_exp.groupby("time_days", as_index=False).agg({
		"density": ["mean", "std", "count"]
	}).reset_index(drop=True)
	df_exp_stats.columns = ["time_days", "density", "std_count", "n_measurements"]

	times_exp = df_exp_stats["time_days"].values
	densities = df_exp_stats["density"].values
	std_counts = df_exp_stats["std_count"].fillna(0).values
	n_measurements = df_exp_stats["n_measurements"].fillna(0).values
	sd_counts = std_counts

	# ========================================================================
	# Scale Data
	# ========================================================================
	# sem_densities = sem_counts * DENSITY_SCALE
	sd_densities = sd_counts * DENSITY_SCALE
	mean_sim_scaled = mean_sim * DENSITY_SCALE
	std_sim_scaled = std_sim * DENSITY_SCALE

	# ========================================================================
	# Create Visualization
	# ========================================================================
	fig, ax = plt.subplots(figsize=(8,8))
	sim_label = r"$\bar{a}_{\mathrm{sim}} \pm \mathrm{SD}$"
	exp_label = r"$\bar{a}_{\mathrm{exp}} \pm \mathrm{SD}$"
	ax.fill_between( time_sim, mean_sim_scaled - std_sim_scaled, mean_sim_scaled + std_sim_scaled,
				 color=SIM_COLOR, alpha=0.25, zorder=2, label=None)
	ax.plot(time_sim, mean_sim_scaled, label=sim_label, color=SIM_COLOR,linewidth=2.5, zorder=3)
	ax.errorbar(times_exp, densities * DENSITY_SCALE, yerr=sd_densities, fmt="o",
			 color=EXP_COLOR, ecolor=EXP_COLOR, capsize=5, markersize=9,
			 label=exp_label, linewidth=2, zorder=4)

	# ========================================================================
	# Configure Axes and Formatting
	# ========================================================================
	ax.set_xlabel("days (post-injury)", fontsize=16)
	ax.set_ylabel(r"Astrocytes' density $[1/{\mu\mathrm{m}}^2]$", fontsize=16)
	ax.set_ylim(-0.2, 7.2)
	ax.text(-0.01, 1.02, r"(10$^{-4}$)", transform=ax.transAxes, fontsize=16)
	ax.minorticks_on()
	ax.tick_params(axis="x", which="major", direction="out", bottom=True, top=False, length=7, labelsize=16)
	ax.tick_params(axis="y", which="major", direction="out", left=True, right=False, length=7, labelsize=16)
	ax.tick_params(axis="x", which="minor", direction="out", bottom=True, top=False, length=4)
	ax.tick_params(axis="y", which="minor", direction="out", left=True, right=False, length=4)
	ax.legend(fontsize=16, loc="best")
	fig.subplots_adjust(left=0.20, right=0.97, bottom=0.14, top=0.97)
	fig.savefig(os.path.join(base, "calibration_astrocytes.png"), dpi=1500)
	plt.close(fig)


def plot_pressure(base):
	sim_file = os.path.join(base, "nt11_per_roi_results.csv")
	exp_file = os.path.join(base, "experimental_pressure_mild.csv")

	# ========================================================================
	# Load Simulation Data
	# ========================================================================
	df_sim = pd.read_csv(sim_file)
	time_sim = df_sim["time"].values

	roi_columns = [col for col in df_sim.columns if col.startswith("ROI") and col.endswith("_avg")]
	roi_data = df_sim[roi_columns].values

	mean_sim = np.mean(roi_data, axis=1)
	std_sim = np.std(roi_data, axis=1)

	# ========================================================================
	# Load Experimental Data
	# ========================================================================
	df_exp = pd.read_csv(exp_file)
	times_exp = df_exp["time_days"].values
	pressure_exp = df_exp["pressure"].values
	std_exp = df_exp["SD"].fillna(0).values

	# ========================================================================
	# Create Visualization
	# ========================================================================
	fig, ax = plt.subplots(figsize=(8,8))
	sim_label = r"$\bar{p}_{\mathrm{sim}} \pm \mathrm{SD}$"
	exp_label = r"$\bar{p}_{\mathrm{exp}} \pm \mathrm{SD}$"
	ax.fill_between(time_sim, mean_sim - std_sim, mean_sim + std_sim, 
				 color=SIM_COLOR, alpha=0.25, zorder=2, label=None)
	ax.plot(time_sim, mean_sim, label=sim_label, color=SIM_COLOR, linewidth=2.5, zorder=3)
	ax.errorbar(times_exp, pressure_exp, yerr=std_exp, fmt="o", 
			 color=EXP_COLOR, ecolor=EXP_COLOR, capsize=5, markersize=9,
			 label=exp_label, linewidth=2, zorder=4)

	# ========================================================================
	# Configure Axes and Formatting
	# ========================================================================
	ax.set_xlabel("days (post-injury)", fontsize=16)
	ax.set_ylabel(r"Pressure $[\mathrm{Pa}]$", fontsize=18)
	ax.set_ylim(-50, 1550)
	ax.set_yticks([0, 300, 600, 900, 1200, 1500])
	ax.minorticks_on()
	ax.tick_params(axis="x", which="major", direction="out", bottom=True, top=False, length=7, labelsize=16)
	ax.tick_params(axis="y", which="major", direction="out", left=True, right=False,length=7, labelsize=16)
	ax.tick_params(axis="x", which="minor", direction="out", bottom=True, top=False,length=4)
	ax.tick_params(axis="y", which="minor", direction="out", left=True, right=False, length=4)
	ax.legend(fontsize=16, loc="best")
	fig.subplots_adjust(left=0.20, right=0.97, bottom=0.14, top=0.97)
	fig.savefig(os.path.join(base, "calibration_pressure.png"), dpi=300)
	plt.close()


def plot_loss_runs(base):
	input_file = Path(base) / "calib_result_de.json"
	output_file = Path(base) / "de_losses_per_run.png"

	with input_file.open("r", encoding="utf-8") as f:
		result = json.load(f)

	losses = np.asarray(result.get("losses", []), dtype=float)
	run_idx = np.arange(1, len(losses) + 1)

	fig, ax = plt.subplots(figsize=(8, 6))
	ax.plot(run_idx, losses, "o-", color=SIM_COLOR, linewidth=2.0, markersize=1)

	ax.set_xlabel("Run", fontsize=16)
	ax.set_ylabel(r"$f(\mathbf{\theta})$", fontsize=16)
	ax.tick_params(axis="both", which="major", labelsize=16)

	fig.tight_layout()
	fig.savefig(output_file, dpi=180)
	plt.close(fig)

# ============================================================================
# DATA LOADING AND PROCESSING
# ============================================================================

def main():
	base = os.path.dirname(os.path.abspath(__file__))
	plot_density(base)
	plot_pressure(base)
	plot_loss_runs(base)


if __name__ == '__main__':
	main()

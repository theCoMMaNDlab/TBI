# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ------------------------------------------------------------------
# Configurations
# ------------------------------------------------------------------

CUSTOM_PURD = [
    "#fef6f7", "#fde9ed", "#fcdbde", "#faccce", "#f7b9bf",
    "#f4a4b0", "#f08fa2", "#ec7894", "#e86187", "#e3497b",
    "#de3070", "#d71965", "#ce035b", "#b80353", "#a1034a",
    "#8a0341", "#720338", "#5b022e", "#440224", "#2d011a",
]
HEATMAP_CMAP = colors.ListedColormap(CUSTOM_PURD, name="CustomPuRd")

TARGET_VALUES = [400, 450, 500]
LINE_COLORS = ["#f4a4b0", "#d71965", "#440224"]

TIME_PLOT_DATA = "Kh_Output.csv"
HEATMAP_DATA = "parametric_results_biot_Kh.csv"

X_COL = "biot"
Y_COL = "Kh"
BASELINE_X = 1.0
BASELINE_Y = 400 

# ------------------------------------------------------------------
# Functions for generating time-series plots and heatmaps
# ------------------------------------------------------------------
def generate_time_series_plots(csv_file):
    df = pd.read_csv(csv_file)
    required_columns = ["Time"]
    for value in TARGET_VALUES:
        required_columns += [
            f"biot=1.0-Kh={value}_NT11",
            f"biot=1.0-Kh={value}_U",
        ]

    plot_specs = [("NT11", "Pressure [Pa]", "time_pressure.png", True),
                  ("U", "Displacement [$\mu$m]", "time_displacement.png", False)]

    for metric_suffix, y_label, output_name, is_pressure_plot in plot_specs:
        plt.figure(figsize=(6, 6))
        for index, value in enumerate(TARGET_VALUES):
            plt.plot(df["Time"], df[f"biot=1.0-Kh={value}_{metric_suffix}"],
                     color=LINE_COLORS[index], linewidth=2, label=rf"$K_{{h}} = {value}$")

        plt.xlabel("days (post-injury)", fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.legend(loc="best", fontsize=14)
        plt.xlim(-0.5, 7.5)

        if is_pressure_plot:
            plt.ylim(-50, 1650)
            plt.yticks([0, 400, 800, 1200, 1600], ["0.0", "400", "800", "1200", "1600"])

        plt.minorticks_on()
        plt.tick_params(axis="both", which="major", labelsize=14)
        plt.tick_params(axis="both", which="minor", labelsize=12)
        plt.tight_layout()
        plt.savefig(f"{output_name}", dpi=500, bbox_inches="tight")
        plt.close()

def generate_heatmaps(csv_file, x_col, y_col, baseline_x_value, baseline_y_value):
    df = pd.read_csv(csv_file)
    missing_axis_columns = [col for col in [x_col, y_col] if col not in df.columns]
    if missing_axis_columns:
        raise ValueError(f"Missing axis columns in {csv_file}: {missing_axis_columns}")

    baseline_row = df.loc[
        (df[x_col] == baseline_x_value) & (df[y_col] == baseline_y_value),
        ["P", "U"]]

    baseline_values = baseline_row.iloc[0]
    for column in ["P", "U"]:
        df[column] = (df[column] - baseline_values[column]).abs()

    heatmaps = [
        ("P", "root mean square error [Pa]", "heatmap_P.png"),
        ("U", r"root mean square error [$\mu$m]", "heatmap_U.png"),
    ]

    for column, label, output_name in heatmaps:
        heat_data = df.pivot(index=y_col, columns=x_col, values=column)
        heat_data = heat_data.sort_index(ascending=False)
        _, ax = plt.subplots(figsize=(6, 6))
        image = ax.imshow(heat_data.values, cmap=HEATMAP_CMAP, aspect="equal", interpolation="nearest")

        ax.set_xticks(range(len(heat_data.columns)))
        ax.set_yticks(range(len(heat_data.index)))
        ax.set_xticklabels(heat_data.columns)
        ax.set_yticklabels(heat_data.index)
        ax.set_ylabel("hydraulic conductivity", fontsize=14)
        ax.set_xlabel("Biot coefficient", fontsize=14)
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center", fontsize=12)
        ax.set_yticklabels([int(float(tick.get_text())) for tick in ax.get_yticklabels()],
                        rotation=0, fontsize=12)

        ax.set_xticks(np.arange(-0.5, len(heat_data.columns), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(heat_data.index), 1), minor=True)
        ax.grid(which="minor", color="white", linestyle="-", linewidth=1.1)
        ax.tick_params(which="minor", bottom=False, left=False)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color("black")
            spine.set_linewidth(1.2)

        ax.set_facecolor("#f3eef0")
        ax.tick_params(axis="both", which="major", length=0)
        colorbar = plt.colorbar(image, ax=ax, orientation="horizontal",
                                shrink=0.72, pad=0.13, fraction=0.07, aspect=25)
        colorbar.set_label(label, fontsize=12)
        colorbar.ax.tick_params(labelsize=11)

        plt.tight_layout()
        plt.savefig(f"{output_name}", dpi=500, bbox_inches="tight")
        plt.close()

# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    generate_time_series_plots(csv_file = TIME_PLOT_DATA)
    generate_heatmaps(csv_file=HEATMAP_DATA, x_col=X_COL, y_col=Y_COL,
                      baseline_x_value=BASELINE_X, baseline_y_value=BASELINE_Y)


if __name__ == "__main__":
    main()

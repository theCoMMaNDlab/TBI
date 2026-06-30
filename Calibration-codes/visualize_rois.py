# Author: Fatemeh Jafarabadi, Apr 2026
# Units: Length [micrometers], Time [days], Stress [Pa]
# ============================================================================
import json
import numpy as np
import matplotlib.pyplot as plt

# Sampling and plotting constants used for readability and reproducibility.
ROI_MAX_ATTEMPTS_FACTOR = 100
SECTOR_ARC_SAMPLES = 200
DOMAIN_ARC_SAMPLES = 400
ROI_WIDTH = 1330
ROI_HEIGHT = 1330
N_ROIS = 3
AXIS_PADDING = 5
ROI_PLOT_DPI = 300
DOMAIN_RADIUS = 5000
INJURY_CENTER_X = 2000
INJURY_CENTER_Y = 4000
INNER_RADIUS = 0
OUTER_RADIUS = 3000
ANGLE_START = 0
ANGLE_END = 360

# ============================================================================
# Injury location definition
# ============================================================================
def is_point_in_injury_region(x, y, injury_center_x, injury_center_y, inner_radius, outer_radius,
                              angle_start, angle_end, domain_radius):
    
    """Return whether a point lies inside the injury-sector region and domain."""
    if y < 0 or x**2 + y**2 > domain_radius**2:
        return False
    
    # Compute distance and angle from injury center
    dx = x - injury_center_x
    dy = y - injury_center_y
    distance = np.sqrt(dx**2 + dy**2)
    angle = np.degrees(np.arctan2(dy, dx))
    
    # Normalize angle to [0, 360)
    angle = angle % 360
    
    # Check radial constraint
    if distance < inner_radius or distance > outer_radius:
        return False
    
    # Check angular constraint
    if angle_start <= angle_end:
        return angle_start <= angle <= angle_end
    else:
        return angle >= angle_start or angle <= angle_end

# ============================================================================
# ROI generation logic
# ============================================================================
def generate_random_rois(injury_center_x=INJURY_CENTER_X, injury_center_y=INJURY_CENTER_Y, inner_radius=INNER_RADIUS,
                         outer_radius=OUTER_RADIUS, angle_start=ANGLE_START, angle_end=ANGLE_END, roi_width=ROI_WIDTH,
                         roi_height=ROI_HEIGHT, n_rois=N_ROIS, out_file="rois_random_10X.json", radius=DOMAIN_RADIUS):
    
    """Generate non-overlapping ROIs inside the injury sector and save to JSON."""
    rois = []
    # Minimum distance between ROI centers (80% of minimum ROI dimension)
    min_dist = min(roi_width, roi_height) * 0.8
    min_dist_sq = min_dist * min_dist
    
    # Maximum attempts to avoid infinite loop
    max_attempts = n_rois * ROI_MAX_ATTEMPTS_FACTOR
    attempts = 0
    
    # Generate ROIs until we have enough non-overlapping ones
    while len(rois) < n_rois and attempts < max_attempts:
        attempts += 1
        # Sample a random point within the injury sector
        r = np.random.uniform(inner_radius, outer_radius)
        angle_rad = np.radians(np.random.uniform(angle_start, angle_end))
        
        # Convert polar to Cartesian coordinates centered at injury location
        cx = injury_center_x + r * np.cos(angle_rad)
        cy = injury_center_y + r * np.sin(angle_rad)
        
        # Get ROI corners for domain boundary check
        half_w = roi_width / 2
        half_h = roi_height / 2
        corners = [
            (cx - half_w, cy - half_h),
            (cx + half_w, cy - half_h),
            (cx - half_w, cy + half_h),
            (cx + half_w, cy + half_h),
        ]
        
        # Check if all corners are within the half-disk domain (x² + y² ≤ R² and y ≥ 0)
        if not all((x**2 + y**2 <= radius**2 and y >= 0) for x, y in corners):
            continue
        
        # Check if ROI center is in injury region (including domain constraint)
        if not is_point_in_injury_region(cx, cy, injury_center_x, injury_center_y, inner_radius, outer_radius, 
                                         angle_start, angle_end, radius):
            continue
        
        # Check for overlap with existing ROIs
        overlap = False
        for (ocx, ocy, _, _) in rois:
            dist_sq = (cx - ocx)**2 + (cy - ocy)**2
            if dist_sq < min_dist_sq:
                overlap = True
                break
        
        if overlap:
            continue
        
        # Accept this ROI
        rois.append([cx, cy, roi_width, roi_height])
    # Save ROIs to JSON file
    with open(out_file, "w") as f:
        json.dump(rois, f, indent=2)
    
    return rois

# ============================================================================
# Visualization 
# ============================================================================

def plot_injury_circle(ax, injury_center_x, injury_center_y, outer_radius, angle_start, 
                       angle_end, domain_radius=DOMAIN_RADIUS):
    """Plot the visible outer arc of the injury sector within the half-disk domain."""
    # Generate points for the outer arc
    theta_outer = np.linspace(np.radians(angle_start), np.radians(angle_end), SECTOR_ARC_SAMPLES)
    x_outer = injury_center_x + outer_radius * np.cos(theta_outer)
    y_outer = injury_center_y + outer_radius * np.sin(theta_outer)
    
    # Filter outer arc to only points within domain (y >= 0 and x² + y² <= R²).
    valid_mask = (y_outer >= 0) & (x_outer**2 + y_outer**2 <= domain_radius**2)
    if np.any(valid_mask):
        x_outer_valid = x_outer[valid_mask]
        y_outer_valid = y_outer[valid_mask]
        ax.plot(x_outer_valid, y_outer_valid, color="firebrick", linewidth=2.0, linestyle='--')

def plot_rois(ax, rois):
    """Draw ROI rectangles and center labels on the provided axes."""
    for i, (cx, cy, w, h) in enumerate(rois):
        rect = plt.Rectangle((cx - w / 2, cy - h / 2), w, h, edgecolor='black', facecolor='none', lw=1.5, linestyle='-')
        ax.add_patch(rect)
        ax.text(cx, cy, str(i + 1), color='black', fontsize=12, ha='center', va='center')

# ============================================================================
# Main
# ============================================================================
def main():
    """Load saved ROIs and render a diagnostic geometry visualization."""

    # Generate random ROIs within the injury region
    roi_file = "rois_random.json"
    # generate_random_rois(injury_center_x=INJURY_CENTER_X, injury_center_y=INJURY_CENTER_Y, inner_radius=INNER_RADIUS,
                        #  outer_radius=OUTER_RADIUS, angle_start=ANGLE_START, angle_end=ANGLE_END,
                        #  roi_width=ROI_WIDTH, roi_height=ROI_HEIGHT, n_rois=N_ROIS, 
                        #  out_file=roi_file, radius=DOMAIN_RADIUS)

    rois = []
    with open(roi_file, 'r') as f:
        rois = json.load(f)
    
    _, ax = plt.subplots(figsize=(6, 6))
    theta = np.linspace(0, np.pi, DOMAIN_ARC_SAMPLES)
    x = DOMAIN_RADIUS * np.cos(theta)
    y = DOMAIN_RADIUS * np.sin(theta)
    ax.plot(x, y, color='k', linewidth=1.5)
    ax.plot([-DOMAIN_RADIUS, DOMAIN_RADIUS], [0, 0], color='k', linewidth=2.0)
    ax.set_aspect('equal')
    for spine in ax.spines.values():
        spine.set_color('gray')
    ax.spines['bottom'].set_color('black')
    ax.set_xlim(-DOMAIN_RADIUS, DOMAIN_RADIUS)
    ax.set_ylim(0, DOMAIN_RADIUS)
    ax.set_xticks([-DOMAIN_RADIUS, -DOMAIN_RADIUS/2, 0, DOMAIN_RADIUS/2, DOMAIN_RADIUS])
    ax.set_yticks([0, DOMAIN_RADIUS/2, DOMAIN_RADIUS])
    ax.set_xlabel(r'$x$ [$\mu$m]', fontdict={'size': 14})
    ax.set_ylabel(r'$y$ [$\mu$m]', fontdict={'size': 14})
    ax.tick_params(axis='both', which='major', labelsize=12)

    # Plot the injury sector region
    plot_injury_circle(ax, INJURY_CENTER_X, INJURY_CENTER_Y, outer_radius=OUTER_RADIUS,
                       angle_start=ANGLE_START, angle_end=ANGLE_END, domain_radius=DOMAIN_RADIUS)

    plot_rois(ax, rois)
    ax.scatter([INJURY_CENTER_X], [INJURY_CENTER_Y], marker='o', s=100, color='firebrick')
    plt.tight_layout()
    plt.savefig("rois_visualization.png", dpi=ROI_PLOT_DPI)
    plt.show()

if __name__ == "__main__":
    main()

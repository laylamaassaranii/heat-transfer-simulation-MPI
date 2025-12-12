import sys
import os
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt


def plot_heatmap(csv_path: str, out_dir: str = "plots", show: bool = True, dpi: int = 200):
    # Read CSV produced by the C++ code
    data = np.loadtxt(csv_path, delimiter=",")

    # Make output directory
    os.makedirs(out_dir, exist_ok=True)

    # Timestamp for unique filenames
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Build output filename
    base = os.path.splitext(os.path.basename(csv_path))[0]
    out_path = os.path.join(out_dir, f"{base}_{timestamp}.png")

    # Plot
    plt.figure(figsize=(6, 5))
    im = plt.imshow(
        data,
        origin="lower",
        cmap="viridis",
        aspect="equal"
    )
    plt.colorbar(im, label="Temperature")
    plt.xlabel("i (x index)")
    plt.ylabel("j (y index)")
    plt.title(f"Final temperature field\n{csv_path}")
    plt.tight_layout()

    # Save BEFORE show()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"Saved figure to: {out_path}")

    if show:
        plt.show()

    plt.close()  # free memory


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python visualization.py <csv_file> [out_dir]")
        print("Example: python visualization.py heat_mpi_large_np4_final.csv plots")
        sys.exit(1)

    csv_file = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) >= 3 else "plots"
    plot_heatmap(csv_file, out_dir=out_dir)

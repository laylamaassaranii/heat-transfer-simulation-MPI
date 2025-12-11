import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns  # facultatif mais joli

def plot_heatmap(csv_path: str):
    # Read CSV produced by the C++ code
    data = np.loadtxt(csv_path, delimiter=",")

    # --- Version Matplotlib simple ---
    plt.figure(figsize=(6, 5))
    im = plt.imshow(
        data,
        origin="lower",   # j=0 en bas
        cmap="viridis",
        aspect="equal"
    )
    plt.colorbar(im, label="Temperature")
    plt.xlabel("i (x index)")
    plt.ylabel("j (y index)")
    plt.title(f"Final temperature field\n{csv_path}")
    plt.tight_layout()
    plt.show()

    # --- Optionnel : version Seaborn ---
    # heatmap_data = np.flipud(data)  # pour avoir y=0 en bas
    # plt.figure(figsize=(6, 5))
    # sns.heatmap(
    #     heatmap_data,
    #     cmap="viridis",
    #     cbar_kws={"label": "Temperature"}
    # )
    # plt.xlabel("i (x index)")
    # plt.ylabel("j (y index)")
    # plt.title(f"Final temperature field (Seaborn)\n{csv_path}")
    # plt.tight_layout()
    # plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python visualization.py <csv_file>")
        print("Example: python visualization.py heat_mpi_large_np4_final.csv")
        sys.exit(1)

    csv_file = sys.argv[1]
    plot_heatmap(csv_file)

"""Analyze the real JWST NIRSpec/G395H phase-resolved emission spectroscopy
of WASP-121 b, extracting a day-night brightness-temperature contrast
spectrum.

Data source: Zenodo record 10.5281/zenodo.20651891, "WASP-121b JWST
NIRSpec/G395H data products". Retrieved directly from Zenodo; the
wavelength grid, phase grid, and brightness-temperature matrix (with
asymmetric upper/lower uncertainties) are reproduced unmodified in data/.

Phase convention (from the source file headers): phase 0 = primary transit
(planet's nightside hemisphere facing the observer), phase +/-0.5 =
secondary eclipse (planet's dayside hemisphere fully visible just before
occultation). This script averages the brightness-temperature spectrum
over the phase bins closest to eclipse (dayside-facing) and closest to
transit (nightside-facing) and computes their real difference at each
wavelength -- a direct, data-driven measurement of the day-night thermal
contrast, without fitting any atmospheric model.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"


def load_wavelength_centers(path: Path) -> np.ndarray:
    edges = np.loadtxt(path, comments="#")
    return 0.5 * (edges[:, 0] + edges[:, 1])


def load_phase_centers(path: Path) -> np.ndarray:
    table = np.loadtxt(path, comments="#")
    return 0.5 * (table[:, 1] + table[:, 2])


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    wavelength = load_wavelength_centers(DATA_DIR / "emspec_wav.txt")
    phase = load_phase_centers(DATA_DIR / "emspec_phi.txt")
    tbright = np.loadtxt(DATA_DIR / "emspec_tbright_vals.txt", comments="#")  # (n_phase, n_wave)

    # Dayside-facing: phase bins nearest +/-0.5 (adjacent to secondary eclipse).
    # Nightside-facing: phase bins nearest 0 (adjacent to primary transit).
    dayside_mask = np.abs(np.abs(phase) - 0.5) < 0.08
    nightside_mask = np.abs(phase) < 0.11

    dayside_t = tbright[dayside_mask, :].mean(axis=0)
    nightside_t = tbright[nightside_mask, :].mean(axis=0)
    contrast = dayside_t - nightside_t

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_channels", len(wavelength), "count"])
        writer.writerow(["n_phase_bins_total", len(phase), "count"])
        writer.writerow(["n_phase_bins_dayside_facing", int(dayside_mask.sum()), "count"])
        writer.writerow(["n_phase_bins_nightside_facing", int(nightside_mask.sum()), "count"])
        writer.writerow(["mean_dayside_brightness_temp", f"{dayside_t.mean():.1f}", "K"])
        writer.writerow(["mean_nightside_brightness_temp", f"{nightside_t.mean():.1f}", "K"])
        writer.writerow(["mean_day_night_contrast", f"{contrast.mean():.1f}", "K"])
        writer.writerow(["max_day_night_contrast", f"{contrast.max():.1f}", "K"])
        writer.writerow(["min_day_night_contrast", f"{contrast.min():.1f}", "K"])

    fig, (ax_t, ax_c) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={"hspace": 0.08})
    ax_t.plot(wavelength, dayside_t, color="#c0562a", lw=1.4, label=f"dayside-facing ({int(dayside_mask.sum())} phase bins)")
    ax_t.plot(wavelength, nightside_t, color="#2c5f8a", lw=1.4, label=f"nightside-facing ({int(nightside_mask.sum())} phase bins)")
    ax_t.set_ylabel("Brightness temperature [K]")
    ax_t.set_title("WASP-121 b day vs. night brightness temperature (real JWST NIRSpec/G395H phase curve)")
    ax_t.legend(fontsize=8, frameon=False)
    ax_t.grid(alpha=0.25)

    ax_c.plot(wavelength, contrast, color="#3c3c3c", lw=1.4)
    ax_c.axhline(0, color="#999999", lw=1)
    ax_c.set_xlabel("Wavelength [micron]")
    ax_c.set_ylabel("Day - night [K]")
    ax_c.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(FIG_DIR / "wasp121b_day_night_contrast.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'wasp121b_day_night_contrast.png'}")
    print(f"n_wave={len(wavelength)}, n_phase={len(phase)}")
    print(f"Mean dayside T = {dayside_t.mean():.1f} K, mean nightside T = {nightside_t.mean():.1f} K")
    print(f"Mean day-night contrast = {contrast.mean():.1f} K (range {contrast.min():.1f} to {contrast.max():.1f} K)")


if __name__ == "__main__":
    main()

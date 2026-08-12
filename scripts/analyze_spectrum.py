"""Analyze the real JWST NIRSpec/G395H phase-resolved emission spectroscopy
of WASP-121 b, extracting a day-night brightness-temperature contrast
spectrum.

Data source: Zenodo record 10.5281/zenodo.20651891, "WASP-121b JWST
NIRSpec/G395H data products". Retrieved directly from Zenodo; the
wavelength grid, phase grid, and brightness-temperature matrix (with
asymmetric upper/lower posterior uncertainties) are reproduced unmodified
in data/.

Phase convention (from the source file headers): phase 0 = primary transit
(planet's nightside hemisphere facing the observer), phase +/-0.5 =
secondary eclipse (planet's dayside hemisphere fully visible just before
occultation). This script combines the brightness-temperature spectrum
over the phase bins closest to eclipse (dayside-facing) and closest to
transit (nightside-facing) with an inverse-variance weighted mean, using
each phase/wavelength point's own posterior uncertainty (averaging the
asymmetric upper/lower bounds into a single sigma), and propagates that
into an uncertainty on the day-night contrast at each wavelength -- rather
than a bare difference of unweighted means with no error bar.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np

plt.style.use(["science", "no-latex"])

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
FIG_DIR = Path(__file__).resolve().parents[1] / "figures"


def load_wavelength_centers(path: Path) -> np.ndarray:
    edges = np.loadtxt(path, comments="#")
    return 0.5 * (edges[:, 0] + edges[:, 1])


def load_phase_centers(path: Path) -> np.ndarray:
    table = np.loadtxt(path, comments="#")
    return 0.5 * (table[:, 1] + table[:, 2])


def weighted_mean_axis0(values: np.ndarray, errors: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Inverse-variance weighted mean down axis 0 (phase), per wavelength
    column. A handful of posterior points in this dataset have zero
    reported uncertainty (a fitting artifact, not a real infinitely
    precise measurement); those get zero weight instead of the
    infinite weight a literal 1/error**2 would otherwise assign them."""
    with np.errstate(divide="ignore"):
        weights = np.where(errors > 0, 1.0 / errors**2, 0.0)
    mean = np.sum(values * weights, axis=0) / np.sum(weights, axis=0)
    mean_error = np.sqrt(1.0 / np.sum(weights, axis=0))
    return mean, mean_error


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)
    wavelength = load_wavelength_centers(DATA_DIR / "emspec_wav.txt")
    phase = load_phase_centers(DATA_DIR / "emspec_phi.txt")
    tbright = np.loadtxt(DATA_DIR / "emspec_tbright_vals.txt", comments="#")  # (n_phase, n_wave)
    tbright_lo = np.loadtxt(DATA_DIR / "emspec_tbright_uncs_lower.txt", comments="#")
    tbright_hi = np.loadtxt(DATA_DIR / "emspec_tbright_uncs_upper.txt", comments="#")
    tbright_err = 0.5 * (tbright_lo + tbright_hi)  # symmetrized posterior sigma per point

    # Dayside-facing: phase bins nearest +/-0.5 (adjacent to secondary eclipse).
    # Nightside-facing: phase bins nearest 0 (adjacent to primary transit).
    dayside_mask = np.abs(np.abs(phase) - 0.5) < 0.08
    nightside_mask = np.abs(phase) < 0.11

    dayside_t, dayside_t_err = weighted_mean_axis0(tbright[dayside_mask, :], tbright_err[dayside_mask, :])
    nightside_t, nightside_t_err = weighted_mean_axis0(tbright[nightside_mask, :], tbright_err[nightside_mask, :])
    contrast = dayside_t - nightside_t
    contrast_err = np.sqrt(dayside_t_err**2 + nightside_t_err**2)

    mean_contrast = np.average(contrast, weights=1.0 / contrast_err**2)
    mean_contrast_err = np.sqrt(1.0 / np.sum(1.0 / contrast_err**2))

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "value", "unit"])
        writer.writerow(["n_wavelength_channels", len(wavelength), "count"])
        writer.writerow(["n_phase_bins_total", len(phase), "count"])
        writer.writerow(["n_phase_bins_dayside_facing", int(dayside_mask.sum()), "count"])
        writer.writerow(["n_phase_bins_nightside_facing", int(nightside_mask.sum()), "count"])
        writer.writerow(["weighted_mean_dayside_brightness_temp", f"{np.average(dayside_t, weights=1/dayside_t_err**2):.1f} +/- {np.sqrt(1/np.sum(1/dayside_t_err**2)):.1f}", "K"])
        writer.writerow(["weighted_mean_nightside_brightness_temp", f"{np.average(nightside_t, weights=1/nightside_t_err**2):.1f} +/- {np.sqrt(1/np.sum(1/nightside_t_err**2)):.1f}", "K"])
        writer.writerow(["weighted_mean_day_night_contrast", f"{mean_contrast:.1f} +/- {mean_contrast_err:.1f}", "K"])
        writer.writerow(["max_day_night_contrast", f"{contrast.max():.1f} +/- {contrast_err[np.argmax(contrast)]:.1f}", "K"])
        writer.writerow(["min_day_night_contrast", f"{contrast.min():.1f} +/- {contrast_err[np.argmin(contrast)]:.1f}", "K"])

    fig, (ax_t, ax_c) = plt.subplots(2, 1, figsize=(9, 7), sharex=True, gridspec_kw={"hspace": 0.08})
    ax_t.fill_between(wavelength, dayside_t - dayside_t_err, dayside_t + dayside_t_err, color="#c0562a", alpha=0.2, lw=0)
    ax_t.plot(wavelength, dayside_t, color="#c0562a", lw=1.4, label=f"dayside-facing ({int(dayside_mask.sum())} phase bins)")
    ax_t.fill_between(wavelength, nightside_t - nightside_t_err, nightside_t + nightside_t_err, color="#2c5f8a", alpha=0.2, lw=0)
    ax_t.plot(wavelength, nightside_t, color="#2c5f8a", lw=1.4, label=f"nightside-facing ({int(nightside_mask.sum())} phase bins)")
    ax_t.set_ylabel("Brightness temperature [K]")
    ax_t.set_title("WASP-121 b day vs. night brightness temperature (JWST NIRSpec/G395H phase curve)")
    ax_t.legend(fontsize=8, frameon=False)
    ax_t.grid(alpha=0.25)

    ax_c.fill_between(wavelength, contrast - contrast_err, contrast + contrast_err, color="#3c3c3c", alpha=0.2, lw=0)
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
    print(f"Weighted mean day-night contrast = {mean_contrast:.1f} +/- {mean_contrast_err:.1f} K (range {contrast.min():.1f} to {contrast.max():.1f} K)")


if __name__ == "__main__":
    main()

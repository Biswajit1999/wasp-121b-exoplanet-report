"""Executable checks on the per-wavelength-column weighted mean and a
regression guard that the pipeline still reproduces the documented
headline numbers when run on the real downloaded data."""

import csv

import numpy as np
import analyze_spectrum as spec


def test_weighted_mean_axis0_matches_hand_computed_case():
    # Two wavelength columns, two phase rows each.
    values = np.array([[1.0, 10.0], [2.0, 20.0]])
    errors = np.array([[1.0, 2.0], [0.5, 1.0]])  # weights [1,4] and [0.25,1]
    mean, mean_err = spec.weighted_mean_axis0(values, errors)
    # Column 0: weighted mean of (1, 2) with weights (1, 4) = 9/5 = 1.8
    assert np.isclose(mean[0], 1.8, rtol=1e-10)
    # Column 1: weighted mean of (10, 20) with weights (0.25, 1) = 22.5/1.25 = 18.0
    assert np.isclose(mean[1], 18.0, rtol=1e-10)


def test_weighted_mean_axis0_zero_error_gets_zero_weight():
    # A point with zero reported uncertainty must not blow up to
    # infinite weight and dominate the mean.
    values = np.array([[5.0], [100.0]])
    errors = np.array([[0.0], [1.0]])
    mean, _ = spec.weighted_mean_axis0(values, errors)
    assert np.isclose(mean[0], 100.0, rtol=1e-6)


def test_pipeline_reproduces_documented_headline_numbers():
    spec.FIG_DIR.mkdir(exist_ok=True)
    spec.main()
    rows = {}
    units = {}
    with (spec.FIG_DIR / "summary_statistics.csv").open() as f:
        for row in csv.DictReader(f):
            rows[row["quantity"]] = row["value"]
            units[row["quantity"]] = row["unit"]

    assert int(rows["n_wavelength_channels"]) == 349
    # The scalar temperature must still be explicitly labeled as
    # non-bolometric wherever it's reported, not silently reverted.
    assert "NOT a bolometric" in units["wavelength_averaged_monochromatic_dayside_Tb_this_script"]

    dayside_val = float(rows["wavelength_averaged_monochromatic_dayside_Tb_this_script"].split(" +/- ")[0])
    contrast_val = float(rows["wavelength_averaged_day_night_contrast_this_script"].split(" +/- ")[0])
    assert abs(dayside_val - 2722.1) < 1.0
    assert abs(contrast_val - 1493.2) < 1.0
    # Published comparison values must be recorded, not silently dropped.
    assert float(rows["paper_bolometric_Tday"].split(" +/- ")[0]) == 2717
    assert float(rows["paper_nightside_Tb_NRS1"].split(" +/- ")[0]) == 926

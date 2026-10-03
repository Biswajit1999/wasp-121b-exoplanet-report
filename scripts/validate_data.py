"""Fail closed if the pinned Zenodo v8 phase-curve products drift."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
EXPECTED_SHA256 = {
    "emspec_phi.txt": "efdb9bb5ab1fc230b216e851d0cf83e1351f90888b969acace505faafb2cc2a2",
    "emspec_wav.txt": "6c4a467dbb36e470b0cf66b377074ffe20e6161d98586ee7899b49d64a5ce869",
    "emspec_tbright_vals.txt": "43de958f6d2377a980f0658354b54cc10154618464a7dd5f42eb8520248a7083",
    "emspec_tbright_uncs_lower.txt": "25b8826c5f462c62334d54734ea8cef56a4abed8b09a315dd6b9fc61ddea4180",
    "emspec_tbright_uncs_upper.txt": "627e1216e28cd9180aeea825f8d3d3840feee7900ca80346a6d3e2812cc8cf41",
    "emspec_fpfs_vals.txt": "e0b233d9510fc92cd101fca3544a7e5eacd3f2fe42171808a586e77dc0c1dc48",
    "emspec_fpfs_uncs.txt": "c1e671c8d3566a543f61117e6b4c2b6d53cf55014e167bb976a46a40981bd354",
}


def canonical_sha256(path: Path) -> str:
    text = path.read_text(encoding="utf-8").replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> None:
    for filename, expected in EXPECTED_SHA256.items():
        actual = canonical_sha256(DATA_DIR / filename)
        if actual != expected:
            raise SystemExit(f"SHA-256 mismatch for {filename}: {actual}")
    phase = np.loadtxt(DATA_DIR / "emspec_phi.txt", comments="#")
    wavelength = np.loadtxt(DATA_DIR / "emspec_wav.txt", comments="#")
    expected_shape = (phase.shape[0], wavelength.shape[0])
    for filename in EXPECTED_SHA256:
        if filename in {"emspec_phi.txt", "emspec_wav.txt"}:
            continue
        array = np.loadtxt(DATA_DIR / filename, comments="#")
        if array.shape != expected_shape or not np.all(np.isfinite(array)):
            raise SystemExit(f"invalid {filename}: shape={array.shape}, expected={expected_shape}")
    print(f"Validated {len(EXPECTED_SHA256)} Zenodo v8 products with shape {expected_shape}")


if __name__ == "__main__":
    main()

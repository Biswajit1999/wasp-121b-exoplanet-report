# Data source

All seven files are downloaded, numerically unmodified, from Zenodo record
**10.5281/zenodo.20651891** ("WASP-121b JWST NIRSpec/G395H data
products"), the data release for Evans-Soma, Sing et al. (2025), "SiO
and a super-stellar C/O ratio in the atmosphere of the giant exoplanet
WASP-121b," *Nature Astronomy*, 9(6), pp.845-861 (arXiv:2506.01771).

- `emspec_wav.txt` — wavelength channel edges (micron)
- `emspec_phi.txt` — orbital phase bin edges (0 = transit, +/-0.5 = eclipse)
- `emspec_tbright_vals.txt` — brightness temperature (K), phase x wavelength matrix
- `emspec_tbright_uncs_upper.txt` / `emspec_tbright_uncs_lower.txt` — asymmetric
  1-sigma uncertainties on the brightness temperature
- `emspec_fpfs_vals.txt` / `emspec_fpfs_uncs.txt` — planet/star flux ratio
  (ppm) and its reported uncertainty, used for detector-band phase curves

The five original files were retrieved on 2026-08-11. The two flux-ratio files
were retrieved on 2026-10-03 through the same API. The record was v8, published
2026-06-12. Numerical comparison confirmed the original five committed arrays
are identical to v8; their raw MD5 differs only because Git converted line
endings. Official v8 MD5 values include `58b0c3c1ffeac3128af0a924f91ab6da`
for `emspec_fpfs_vals.txt` and `a91c1e481b1da24939cf668d9e3100ac`
for `emspec_fpfs_uncs.txt`.

Each file's own header comments (lines starting with `#`) describe its
exact format; see the files themselves for full detail.

`scripts/validate_data.py` pins LF-normalized SHA-256 identities for all seven
files and validates their expected 36 phase-bin × 349 wavelength-channel
matrix shape. This normalization makes validation stable across Git checkouts
on Windows and Linux.

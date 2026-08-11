# Data source

All five files are downloaded, unmodified, from Zenodo record
**10.5281/zenodo.20651891** ("WASP-121b JWST NIRSpec/G395H data products"):

- `emspec_wav.txt` — wavelength channel edges (micron)
- `emspec_phi.txt` — orbital phase bin edges (0 = transit, +/-0.5 = eclipse)
- `emspec_tbright_vals.txt` — brightness temperature (K), phase x wavelength matrix
- `emspec_tbright_uncs_upper.txt` / `emspec_tbright_uncs_lower.txt` — asymmetric
  1-sigma uncertainties on the brightness temperature

Retrieved: 2026-08-11, via `https://zenodo.org/api/records/20651891`.

Each file's own header comments (lines starting with `#`) describe its
exact format; see the files themselves for full detail.

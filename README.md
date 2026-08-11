# WASP-121 b — Exoplanet Atmosphere Report

An ultra-hot Jupiter on a 1.27-day orbit, tidally distorted by its host
star, with a real measured day-night temperature contrast exceeding 1500 K.
This repo extracts that contrast directly from a real JWST phase-resolved
emission spectroscopy dataset.

**[Open the full report](index.html)** (open locally in a browser, or serve
with `python -m http.server` from this directory).

## What's real here

- **System parameters** — queried live from the NASA Exoplanet Archive TAP
  service (`pscomppars` table).
- **Phase-resolved emission spectrum** — real reduced JWST NIRSpec/G395H
  data: brightness temperature as a function of both wavelength (349
  channels, 2.7-5.2 microns) and orbital phase (36 bins spanning almost a
  full orbit), released publicly on Zenodo
  ([10.5281/zenodo.20651891](https://doi.org/10.5281/zenodo.20651891)).
- **Analysis** — `scripts/analyze_spectrum.py` averages the real
  brightness-temperature spectrum over phase bins nearest secondary eclipse
  (dayside-facing) and nearest primary transit (nightside-facing), and
  computes their real difference at each wavelength — no atmospheric model
  fitting involved, just a direct comparison of the observed data. Run it
  yourself:

  ```bash
  pip install -r requirements.txt
  python scripts/analyze_spectrum.py
  ```

## Repository structure

```text
index.html              the report webpage
data/                    real JWST NIRSpec/G395H phase-curve data (Zenodo)
scripts/analyze_spectrum.py   real day/night phase-averaging analysis
figures/                 generated plot + summary_statistics.csv
```

## Key finding this repo shows directly

Mean dayside brightness temperature 2737 K vs. mean nightside 1176 K — a
real, measured day-night contrast of about 1561 K (range 1252-2195 K across
wavelength). This is a direct observational signature of very inefficient
heat redistribution, consistent with the extreme irradiation this planet
receives on its 1.27-day orbit.

## References

1. Delrez, L. et al., 2016. WASP-121 b: a hot Jupiter close to tidal
   disruption transiting an active F star. *Monthly Notices of the Royal
   Astronomical Society*, 458(4), pp.4025-4043.
2. Evans, T.M. et al., 2017. An ultrahot gas-giant exoplanet with a
   stratosphere. *Nature*, 548, pp.58-61.
3. Evans, T.M. et al., 2018. Detection of H2O and Evidence for TiO/VO in an
   Ultra-Hot Exoplanet Atmosphere. *The Astrophysical Journal Letters*, 822,
   L4.
4. Zenodo record
   [10.5281/zenodo.20651891](https://doi.org/10.5281/zenodo.20651891),
   "WASP-121b JWST NIRSpec/G395H data products."
5. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)

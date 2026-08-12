# WASP-121 b — Exoplanet Atmosphere Report

An ultra-hot Jupiter on a 1.27-day orbit, tidally distorted by its host
star, with a measured day-night temperature contrast of roughly 1490 K.
This repo extracts that contrast from a JWST phase-resolved emission
spectroscopy dataset, propagating each point's own posterior
uncertainty rather than reporting a bare difference of averages.

**[Open the full report](index.html)** (open locally in a browser, or serve
with `python -m http.server` from this directory).

## Data sources

- **System parameters** — from the NASA Exoplanet Archive TAP
  service (`pscomppars` table).
- **Phase-resolved emission spectrum** — reduced JWST NIRSpec/G395H
  data: brightness temperature (with asymmetric posterior uncertainty)
  as a function of both wavelength (349 channels, 2.7-5.2 microns) and
  orbital phase (36 bins spanning almost a full orbit), from
  Evans-Soma, Sing et al. (2025), released publicly on Zenodo
  ([10.5281/zenodo.20651891](https://doi.org/10.5281/zenodo.20651891)).
- **Analysis** — `scripts/analyze_spectrum.py` combines the brightness-
  temperature spectrum over phase bins nearest secondary eclipse
  (dayside-facing) and nearest primary transit (nightside-facing) with
  an inverse-variance weighted mean, using each point's own symmetrized
  posterior uncertainty, and propagates that into an error bar on the
  day-night contrast at every wavelength. Run it yourself:

  ```bash
  pip install -r requirements.txt
  python scripts/analyze_spectrum.py
  ```

## Repository structure

```text
index.html              the report webpage
data/                    JWST NIRSpec/G395H phase-curve data (Zenodo)
scripts/analyze_spectrum.py   day/night phase-averaging with propagated uncertainty
figures/                 generated plot + summary_statistics.csv
```

## What the numbers show

Weighted mean dayside brightness temperature 2751 ± 3 K vs. nightside
1252 ± 2 K — a day-night contrast of about 1493 ± 4 K (range roughly
1235-1804 K across wavelength). This is a direct observational
signature of very inefficient heat redistribution, consistent with the
extreme irradiation this planet receives on its 1.27-day orbit. It's a
data comparison, not a model fit, and doesn't by itself constrain wind
speeds or redistribution efficiency the way a general-circulation-model
comparison would.

## Limitations

One posterior point in the released dataset has zero reported
uncertainty — almost certainly a fitting artifact rather than an
infinitely precise measurement. The script assigns it zero weight
instead of letting a literal 1/error² blow up into an infinite weight
at that point.

## References

1. Delrez, L. et al., 2016. WASP-121 b: a hot Jupiter close to tidal
   disruption transiting an active F star. *Monthly Notices of the Royal
   Astronomical Society*, 458(4), pp.4025-4043.
2. Evans, T.M. et al., 2017. An ultrahot gas-giant exoplanet with a
   stratosphere. *Nature*, 548, pp.58-61.
3. Evans, T.M. et al., 2018. Detection of H2O and Evidence for TiO/VO in an
   Ultra-Hot Exoplanet Atmosphere. *The Astrophysical Journal Letters*, 822,
   L4.
4. Evans-Soma, T.M., Sing, D.K. et al., 2025. SiO and a super-stellar
   C/O ratio in the atmosphere of the giant exoplanet WASP-121b.
   *Nature Astronomy*, 9(6), pp.845-861 (arXiv:2506.01771).
5. Zenodo record
   [10.5281/zenodo.20651891](https://doi.org/10.5281/zenodo.20651891),
   "WASP-121b JWST NIRSpec/G395H data products."
6. NASA Exoplanet Archive, <https://exoplanetarchive.ipac.caltech.edu/>.

## Author

Biswajit Jana — [Portfolio](https://biswajit1999.github.io/Biswajit_Jana.github.io/) · [GitHub](https://github.com/Biswajit1999) · [LinkedIn](https://www.linkedin.com/in/biswajit-jana-27011a151/) · [ORCID](https://orcid.org/0009-0002-2411-1891)

# Fixture cache

`small-tial` is a reconstruction of the package's supplied fitted modern
ACEpotentials reference (never refitted). Its recipe/resolved export is copied
into this experiment's namespace; original `data/ace_ref` files are unchanged.
`medium-ti` and `large-tial` are deterministic synthetic performance fixtures
created by `generate_modern.jl`, under this repository's code license. Their
nonzero coefficients have no claimed physical accuracy.

`ti-resolved.toml` and `tialv-resolved.toml` are analytic conversions of the
published fitted Ti and April08 Ti–Al–V models by Connor Allen and Albert Bartók,
*Multi-phase dataset for bulk Ti and the Ti-6Al-4V alloy*, Scientific Data 12,
1925 (2025), DOI 10.1038/s41597-025-05302-3. Source archive:
https://doi.org/10.5281/zenodo.15114121 ; license CC BY 4.0:
https://creativecommons.org/licenses/by/4.0/ .

UK Ministry of Defence © British Crown Owned Copyright 2024/AWE-NST

Original archive `README.txt` and `record.json` are retained here. Copyright and
attribution remain with the original authors and applicable British Crown/AWE
notices; these artifacts are separate from the repository code license.
Conversions expand the complex PI coordinates into real products and retain exact
analytic radial recurrences, source coefficients, species, pair/reference terms
and source hashes. They are not refits or splined PACE exports. Regenerate with
`benchmark/gpu_ace/{acquire.py,export_legacy.jl}` in the pinned isolated Julia
1.10.12 / ACEpotentials 0.6.6 environment. JSON source SHA-256 pins are in
`acquire.py`; byte-identical originals remain in `/tmp/pfw-stage1-legacy-inputs`
for this run and are intentionally reacquirable rather than duplicated here.
See `benchmark/gpu_ace/THIRD_PARTY.md` for the original conversion record.

The `*-oracle.toml` files in this directory are **initial exploratory** modern
references collected during setup. Final uncontended oracle measurements live in
`../oracles-final`; benchmark commands explicitly select that directory. Both
sets retain complete numerical references and actual sample counts.

# Benchmark model data

The optional Ti and Ti–Al–V model artifacts are by Connor Allen and Albert
Bartók, associated with *Multi-phase dataset for bulk Ti and the Ti-6Al-4V
alloy*, Scientific Data 12, 1925 (2025), DOI 10.1038/s41597-025-05302-3.
Source archive: https://doi.org/10.5281/zenodo.15114121.

The archive description states: “UK Ministry of Defence © British Crown Owned
Copyright 2024/AWE-NST”. This notice is retained verbatim with the record metadata.

Archive metadata retrieved 2026-11-09 declares CC BY 4.0:
https://creativecommons.org/licenses/by/4.0/. Data copyright and attribution
remain with the original authors, including applicable British Crown/AWE notices;
these artifacts are separate from this repository's code license. Acquisition
retains the original README, complete record metadata and byte-identical JSONs.
No refitting, coefficient replacement or conversion is performed by acquire.py.
The selected alloy fit is April08 (not May30). SHA-256 pins are in acquire.py;
all three acquired files also matched the archive's MD5 metadata during inspection.
No GAP models are downloaded.

The source JSONs store legacy ACE1/JuLIP precontracted complex PI potentials,
not the TiAl recipe/resolved schema. An independent isolated ACEpotentials 0.6.6
process loads both originals successfully. `export_legacy.jl` performs an analytic, non-refitting conversion to the existing
resolved schema with a registered legacy polynomial callable and an algebraic
complex-to-real product expansion. The fitted source PI coefficient order is
retained, with explicit species offsets and directed half-pair weights. Independent
CPU and RX 6800 XT E/F/V comparisons were run in separate environments. Original
invariant fitting coordinates are absent from these precontracted artifacts;
individual PI rows are not presented as invariants. Converted large model files
are generated in the caller-selected cache, not shipped or downloaded at load.
The package's original TiAl reference files remain byte-for-byte untouched.

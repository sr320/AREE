# Metabolomics raw-mode fixture (SIMULATED)

Synthetic, for CI only. Generated with a fixed seed (20261016) to exercise the
targeted-panel options of `workflows/metabolomics` in `raw_reanalysis` mode:

- two platforms (`LC`, `GC`), each with one internal standard
  (`is_internal_standard: true`) and its own per-sample loading, which
  `--metabolomics.normalization internal_standard` must remove;
- 0.0001 below-detection placeholders (`--metabolomics.detection_limit 0.001`):
  `LC_treatonly` is detected only in treatment and `LC_lod` nowhere, so both
  must be reported as not tested, never given a fold change;
- 12 samples in 6 tanks (two families per tank), fitted as 3 vs 3 tank units
  with `--metabolomics.replicate_unit tank`.

Built-in effects (treatment vs control, log2): `LC_up` +2.0, `LC_down` -1.5,
`GC_up` +1.0; the rest are null.

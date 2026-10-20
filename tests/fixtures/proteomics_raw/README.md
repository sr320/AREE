# Proteomics raw-mode fixture (SIMULATED)

Synthetic, for CI only. Generated with a fixed seed (20261017) to exercise
`workflows/proteomics` in `raw_reanalysis` mode on a spectral-count (NSAF-like)
matrix, run with `--proteomics.input_scale linear --proteomics.normalization none
--proteomics.quantification spectral_count`:

- linear values far below 1, which the `auto` scale heuristic would wrongly
  take for log2;
- zeros meaning "not detected": `PX_treatonly` is detected only in treatment
  and `PX_absent` nowhere, so both must be reported as not tested, never given
  a fold change; every tenth null protein has one dropout;
- 10 matrix columns of which the sample sheet names 8 (4 control, 4
  treatment); `other_1` and `other_2` must be left out of the comparison.

Built-in effects (treatment vs control, log2): `PX_up` +2.0, `PX_down` -1.5;
the rest are null.

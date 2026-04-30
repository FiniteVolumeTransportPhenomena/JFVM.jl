# Benchmarks

This folder contains lightweight performance gates used by CI.

- Run manually:

  julia --project benchmarks/perf_gate.jl

The gate checks median time and median allocation for representative 2D kernels:

- linearMean!
- diffusionTerm
- convectionTerm
- transientTerm

Thresholds are intentionally conservative and should be updated when performance
improves materially or when benchmark scope changes.

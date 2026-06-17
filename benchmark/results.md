# Gradient benchmark results

Minimum wall-clock time to evaluate the gradient of each objective once (after warmup),
for the three reverse/forward AD backends and central finite differences, alongside the
primal (function-value) time.  The number in parentheses is the gradient time divided by
the primal time.

These measure *steady-state* (already-compiled) gradient evaluation; the one-time AD
compilation cost is discussed in the "Performance notes" section of the top-level README.

Reproduce with:

```sh
julia --project=benchmark benchmark/benchmarks.jl
```

## Environment

- Julia 1.11.9, 4 × Intel(R) Xeon(R) @ 2.80GHz, Linux x86-64
- DifferentiationInterface v0.7.18, Enzyme v0.13.157, Mooncake v0.5.31,
  FiniteDifferences v0.12.34, StaticArrays v1.9.18


| case | primal | Enzyme reverse | Enzyme forward | Mooncake | FiniteDifferences |
|---|---|---|---|---|---|
| level Ball3 / x | 4.983 ns | 43.600 ns (8.7x) | 208.588 ns (41.9x) | 95.014 ns (19.1x) | 10.227 μs (2052.4x) |
| level Cuboid3 / x | 9.560 ns | 74.954 ns (7.8x) | 220.192 ns (23.0x) | 192.117 ns (20.1x) | 10.073 μs (1053.7x) |
| level Polygon8 / x | 34.698 ns | 107.256 ns (3.1x) | 261.308 ns (7.5x) | 348.190 ns (10.0x) | 7.611 μs (219.3x) |
| level Cylinder / x | 29.956 ns | 72.993 ns (2.4x) | 225.533 ns (7.5x) | 188.598 ns (6.3x) | 10.572 μs (352.9x) |
| surfpt Ball3 / x | 30.027 ns | 59.835 ns (2.0x) | 235.678 ns (7.8x) | 172.052 ns (5.7x) | 10.407 μs (346.6x) |
| surfpt Cuboid3 / x | 34.059 ns | 99.698 ns (2.9x) | 268.272 ns (7.9x) | 713.497 ns (20.9x) | 10.829 μs (317.9x) |
| surfpt Polygon8 / x | 63.033 ns | 107.605 ns (1.7x) | 302.639 ns (4.8x) | 860.695 ns (13.7x) | 9.347 μs (148.3x) |
| surfpt Cylinder / x | 39.161 ns | 120.098 ns (3.1x) | 274.515 ns (7.0x) | 337.577 ns (8.6x) | 11.313 μs (288.9x) |
| volfrac Ball3 / x | 5.465 μs | 16.384 μs (3.0x) | 12.161 μs (2.2x) | 22.526 μs (4.1x) | 249.685 μs (45.7x) |
| volfrac Cuboid3 / x | 5.619 μs | 16.748 μs (3.0x) | 12.301 μs (2.2x) | 26.759 μs (4.8x) | 248.968 μs (44.3x) |
| level Ball3 / params | 5.497 ns | 49.862 ns (9.1x) | 216.187 ns (39.3x) | 94.823 ns (17.2x) | 13.315 μs (2422.2x) |
| level Ellipsoid2 / params | 61.936 ns | 401.260 ns (6.5x) | 721.545 ns (11.6x) | 288.879 ns (4.7x) | 19.168 μs (309.5x) |
| surfpt Polygon5 / params | 5.395 μs | 7.292 μs (1.4x) | 10.487 μs (1.9x) | 19.950 μs (3.7x) | 737.834 μs (136.8x) |
| surfpt Cylinder / params | 159.981 ns | 461.222 ns (2.9x) | 1.392 μs (8.7x) | 1.018 μs (6.4x) | 41.351 μs (258.5x) |

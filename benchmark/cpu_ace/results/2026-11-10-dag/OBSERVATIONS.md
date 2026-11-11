# DAG measured observations

All timings below are warmed medians in milliseconds. Ratios are prefix / DAG. Raw TOML retains every sample, extrema, GC observations and allocation volumes.

| Model | Atoms | Precision | Polynomial prepared | Prefix prepared | DAG prepared | Gain | Prefix public | DAG public | Gain |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|
| small-tial | 128 | Float64 | 1.784 | 1.340 | 1.269 | 1.06× | 3.023 | 2.932 | 1.03× |
| small-tial | 128 | Float32 | 1.665 | 1.235 | 1.164 | 1.06× | 1.520 | 1.448 | 1.05× |
| small-tial | 1024 | Float64 | 14.276 | 10.671 | 10.123 | 1.05× | 25.626 | 24.947 | 1.03× |
| small-tial | 1024 | Float32 | 13.377 | 9.931 | 9.403 | 1.06× | 12.820 | 12.166 | 1.05× |
| medium-ti | 128 | Float64 | 2.331 | 1.644 | 1.540 | 1.07× | 3.144 | 3.052 | 1.03× |
| medium-ti | 128 | Float32 | 2.108 | 1.500 | 1.411 | 1.06× | 1.775 | 1.699 | 1.04× |
| medium-ti | 1024 | Float64 | 18.022 | 13.137 | 12.346 | 1.06× | 26.793 | 25.865 | 1.04× |
| medium-ti | 1024 | Float32 | 16.999 | 12.102 | 11.414 | 1.06× | 14.791 | 13.920 | 1.06× |
| large-tial | 128 | Float64 | 76.847 | 17.821 | 8.680 | 2.05× | 19.513 | 10.333 | 1.89× |
| large-tial | 128 | Float32 | 73.030 | 18.329 | 8.228 | 2.23× | 18.679 | 8.610 | 2.17× |
| large-tial | 1024 | Float64 | 594.411 | 139.909 | 68.453 | 2.04× | 155.992 | 83.257 | 1.87× |
| large-tial | 1024 | Float32 | 572.949 | 138.496 | 64.506 | 2.15× | 141.579 | 68.063 | 2.08× |
| ti | 128 | Float64 | 75.336 | 25.153 | 16.083 | 1.56× | 28.323 | 19.401 | 1.46× |
| ti | 128 | Float32 | 74.076 | 24.285 | 15.062 | 1.61× | 24.607 | 15.515 | 1.59× |
| ti | 1024 | Float64 | 626.939 | 193.468 | 122.529 | 1.58× | 224.799 | 151.518 | 1.48× |
| ti | 1024 | Float32 | 587.077 | 191.120 | 115.272 | 1.66× | 197.060 | 120.629 | 1.63× |
| tialv | 128 | Float64 | 219.731 | 60.771 | 31.312 | 1.94× | 65.070 | 35.320 | 1.84× |
| tialv | 128 | Float32 | 220.971 | 63.017 | 29.396 | 2.14× | 63.768 | 30.394 | 2.10× |
| tialv | 1024 | Float64 | 1794.081 | 524.723 | 243.267 | 2.16× | 556.941 | 275.858 | 2.02× |
| tialv | 1024 | Float32 | 2456.522 | 501.383 | 229.721 | 2.18× | 508.571 | 235.431 | 2.16× |
| small-tial | 11664 | Float64 | — | 121.379 | 115.235 | 1.05× | 313.353 | 313.349 | 1.00× |
| small-tial | 11664 | Float32 | — | 113.434 | 106.800 | 1.06× | 145.568 | 139.497 | 1.04× |
| large-tial | 11664 | Float64 | — | 1663.012 | 775.673 | 2.14× | 1845.275 | 975.139 | 1.89× |
| large-tial | 11664 | Float32 | — | 1590.994 | 737.386 | 2.16× | 1596.433 | 779.988 | 2.05× |

## Fresh calculator plus complete public call

| Model | Atoms | Precision | Prefix ms | DAG ms | Gain |
|---|---:|---|---:|---:|---:|
| small-tial | 128 | Float64 | 3.075 | 3.031 | 1.01× |
| small-tial | 128 | Float32 | 1.568 | 1.507 | 1.04× |
| small-tial | 1024 | Float64 | 25.765 | 25.090 | 1.03× |
| small-tial | 1024 | Float32 | 12.947 | 12.248 | 1.06× |
| medium-ti | 128 | Float64 | 3.217 | 3.110 | 1.03× |
| medium-ti | 128 | Float32 | 1.827 | 1.752 | 1.04× |
| medium-ti | 1024 | Float64 | 26.938 | 26.098 | 1.03× |
| medium-ti | 1024 | Float32 | 15.033 | 14.180 | 1.06× |
| large-tial | 128 | Float64 | 20.914 | 12.227 | 1.71× |
| large-tial | 128 | Float32 | 20.189 | 9.973 | 2.02× |
| large-tial | 1024 | Float64 | 157.932 | 84.434 | 1.87× |
| large-tial | 1024 | Float32 | 143.069 | 69.442 | 2.06× |
| ti | 128 | Float64 | 29.639 | 20.561 | 1.44× |
| ti | 128 | Float32 | 25.935 | 16.916 | 1.53× |
| ti | 1024 | Float64 | 228.310 | 153.732 | 1.49× |
| ti | 1024 | Float32 | 198.965 | 122.348 | 1.63× |
| tialv | 128 | Float64 | 71.974 | 42.250 | 1.70× |
| tialv | 128 | Float32 | 70.326 | 37.363 | 1.88× |
| tialv | 1024 | Float64 | 567.418 | 286.142 | 1.98× |
| tialv | 1024 | Float32 | 495.167 | 241.516 | 2.05× |
| small-tial | 11664 | Float64 | 322.580 | 312.851 | 1.03× |
| small-tial | 11664 | Float32 | 150.402 | 142.064 | 1.06× |
| large-tial | 11664 | Float64 | 1850.498 | 975.889 | 1.90× |
| large-tial | 11664 | Float32 | 1613.589 | 787.655 | 2.05× |

## Product diagnostics at 128 atoms

| Model | Precision | Polynomial product/reverse | Prefix product/reverse | DAG forward | DAG reverse | DAG combined | Product gain | Readout |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| small-tial | Float64 | 0.004205 | 0.000940 | 0.000140 | 0.000270 | 0.000390 | 2.41× | 0.000030 |
| small-tial | Float32 | 0.004180 | 0.000950 | 0.000140 | 0.000260 | 0.000390 | 2.44× | 0.000025 |
| medium-ti | Float64 | 0.006100 | 0.001410 | 0.000190 | 0.000375 | 0.000520 | 2.71× | 0.000030 |
| medium-ti | Float32 | 0.005850 | 0.001320 | 0.000200 | 0.000375 | 0.000520 | 2.54× | 0.000025 |
| large-tial | Float64 | 0.538077 | 0.104471 | 0.009785 | 0.023110 | 0.032845 | 3.18× | 0.000870 |
| large-tial | Float32 | 0.537447 | 0.106950 | 0.009660 | 0.022526 | 0.031711 | 3.37× | 0.000415 |
| ti | Float64 | 0.499727 | 0.105331 | 0.011760 | 0.025851 | 0.037551 | 2.81× | 0.000970 |
| ti | Float32 | 0.498022 | 0.106216 | 0.011600 | 0.025030 | 0.036430 | 2.92× | 0.000530 |
| tialv | Float64 | 1.592035 | 0.349327 | 0.035855 | 0.084726 | 0.115450 | 3.03× | 0.003095 |
| tialv | Float32 | 1.638346 | 0.374036 | 0.034550 | 0.077385 | 0.112556 | 3.32× | 0.001450 |

## Structural work (both precisions share indices)

| Model | Leaves | Original products | Orders | Unique terminals | Multiply nodes | Intermediate only | Depth | Sharing ratio | Prefix fwd/rev multiplies | DAG fwd/rev multiplies |
|---|---:|---:|---|---:|---:|---:|---:|---:|---|---|
| small-tial | 27 | 167 | 1: 12, 2: 63, 3: 92 | 167 | 167 | 12 | 2 | 1.48 | 414 / 828 | 167 / 334 |
| medium-ti | 43 | 230 | 1: 10, 2: 81, 3: 139 | 230 | 257 | 37 | 2 | 1.40 | 589 / 1178 | 257 / 514 |
| large-tial | 165 | 13919 | 1: 24, 2: 630, 3: 4264, 4: 9001 | 13919 | 16229 | 2334 | 3 | 2.23 | 50080 / 100160 | 16229 / 32458 |
| ti | 329 | 16444 | 1: 26, 2: 994, 3: 15424 | 16444 | 19705 | 3287 | 2 | 1.62 | 48286 / 96572 | 19705 / 39410 |
| tialv | 441 | 51379 | 1: 63, 2: 2880, 3: 48436 | 51379 | 57160 | 5844 | 2 | 1.75 | 151131 / 302262 | 57160 / 114320 |

## Construction and memory at 128 atoms, one slot

| Model | Precision | Plan build ms | Plan build allocated bytes | Executor build ms | Plan array payload | DAG scratch | Retained polynomial | Total slot scratch | Fresh prefix public ms | Fresh DAG public ms |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| small-tial | Float64 | 0.0106 | 67640 | 0.0188 | 7320 | 3104 | 256 | 33328 | 3.075 | 3.031 |
| small-tial | Float32 | 0.0112 | 67640 | 0.0170 | 7320 | 1552 | 128 | 16664 | 1.568 | 1.507 |
| medium-ti | Float64 | 0.0147 | 109272 | 0.0244 | 10664 | 4800 | 256 | 52744 | 3.217 | 3.110 |
| medium-ti | Float32 | 0.0148 | 109272 | 0.0236 | 10664 | 2400 | 128 | 26372 | 1.827 | 1.752 |
| large-tial | Float64 | 1.1487 | 7333584 | 1.2732 | 771656 | 262304 | 400 | 484368 | 20.914 | 12.227 |
| large-tial | Float32 | 1.2229 | 7333584 | 1.3805 | 771656 | 131152 | 200 | 242184 | 20.189 | 9.973 |
| ti | Float64 | 1.2946 | 7527664 | 1.3706 | 833120 | 320544 | 256 | 606152 | 29.639 | 20.561 |
| ti | Float32 | 1.1903 | 7527664 | 1.3732 | 833120 | 160272 | 128 | 303076 | 25.935 | 16.916 |
| tialv | Float64 | 4.4594 | 25860176 | 6.6771 | 2534640 | 921616 | 256 | 1516944 | 71.974 | 42.250 |
| tialv | Float32 | 5.4802 | 25860176 | 5.7517 | 2534640 | 460808 | 128 | 758472 | 70.326 | 37.363 |

## Numerical disagreement, max across measured sizes

| Model | Precision | DAG vs polynomial E eV | F eV/Å | W eV | DAG vs independent Float64 E eV | F eV/Å | W eV |
|---|---|---:|---:|---:|---:|---:|---:|
| small-tial | Float64 | 0 | 5.5511e-15 | 1.819e-11 | 1.4901e-08 | 3.8425e-13 | 1.4188e-10 |
| small-tial | Float32 | 0 | 3.0249e-06 | 0.0039062 | 0.12394 | 1.8019e-05 | 0.054709 |
| medium-ti | Float64 | 0 | 3.2196e-15 | 4.5475e-13 | 2.7285e-12 | 2.7145e-13 | 3.7289e-11 |
| medium-ti | Float32 | 0 | 1.6391e-06 | 0.00073242 | 0.00023119 | 4.6495e-05 | 0.00092396 |
| large-tial | Float64 | 0 | 1.2421e-15 | 1.1369e-12 | 5.6843e-14 | 4.2841e-14 | 2.4784e-11 |
| large-tial | Float32 | 0 | 6.1095e-07 | 0.00085449 | 0.00013505 | 1.1672e-05 | 0.0074364 |
| ti | Float64 | 0 | 2.5043e-14 | 1.819e-12 | 1.7462e-10 | 2.634e-14 | 3.6948e-12 |
| ti | Float32 | 0 | 1.277e-05 | 0.0032959 | 0.0080703 | 1.7501e-05 | 0.0016696 |
| tialv | Float64 | 0 | 2.4772e-14 | 9.0949e-13 | 3.7835e-10 | 2.1892e-14 | 2.3874e-12 |
| tialv | Float32 | 0 | 1.3361e-05 | 0.0013428 | 0.19189 | 9.5877e-06 | 0.0012497 |

## Controlled physical-core scaling, 1,024 atoms, Float64

| Cores | Model | Prefix prepared ms | DAG prepared ms | Gain | Prefix public ms | DAG public ms | All-slot scratch bytes | Shared force bytes | DAG prepared host allocation bytes |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | small-tial | 10.989 | 10.205 | 1.08× | 25.923 | 24.994 | 33328 | 24576 | 86864 |
| 1 | large-tial | 139.537 | 68.523 | 2.04× | 155.572 | 83.533 | 484368 | 24576 | 311344 |
| 4 | small-tial | 5.731 | 5.542 | 1.03× | 21.397 | 21.485 | 133312 | 24576 | 236288 |
| 4 | large-tial | 55.215 | 35.911 | 1.54× | 72.855 | 54.145 | 1937472 | 24576 | 460768 |
| 16 | small-tial | 5.958 | 4.565 | 1.31× | 19.087 | 19.100 | 533248 | 24576 | 781456 |
| 16 | large-tial | 38.484 | 21.065 | 1.83× | 58.160 | 37.857 | 7749888 | 24576 | 1005936 |

## Immutable stage-one versus same-checkout prefix

| Model | Atoms | Immutable prepared ms | Integrated prefix ms | Immutable public ms | Integrated prefix public ms | Immutable fresh public ms | Integrated fresh prefix ms |
|---|---:|---:|---:|---:|---:|---:|---:|
| small-tial | 128 | 1.332 | 1.340 | 2.985 | 3.023 | 3.044 | 3.075 |
| small-tial | 1024 | 10.695 | 10.671 | 25.572 | 25.626 | 25.800 | 25.765 |
| large-tial | 128 | 18.208 | 17.821 | 19.944 | 19.513 | 20.091 | 20.914 |
| large-tial | 1024 | 140.958 | 139.909 | 156.013 | 155.992 | 155.925 | 157.932 |

612 trials audited; 0 under-sampled. Counts and full extrema are in summary.json and source TOML.

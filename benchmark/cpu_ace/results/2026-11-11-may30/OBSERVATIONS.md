# Complete May30 observations

All energies and virials are in eV; forces are in eV/Å. Energy error is per atom.

## bcc-128

Composition in Ti/Al/V order: [42, 42, 44]; directed neighbors: 7334. Face heights: [13.73145162018379, 13.53006829469225, 13.681600000000001] Å; minimum distance: 2.8778884752 Å. Verified center-1 V–V nearest image: atom 2, shift [0, 0, 0], distance 2.89633634636 Å.

Cell (column-major): `[13.735999999999999, 0.0, 0.0, 0.31279999999999997, 13.532, 0.0, -0.14959999999999998, 0.23120000000000002, 13.6816]`. Geometry SHA256: `ff7c425703e42ca3fc532116b2738af99efa27450c995e9e11add162472b41bb`.

| Comparison | Energy/atom | Force component | Virial component |
|---|---:|---:|---:|
| ET_conversion | 2.2737e-13 | 1.1324e-14 | 1.1084e-12 |
| dag_prepared | 9.0949e-13 | 2.7898e-12 | 1.2619e-11 |
| dag_public | 9.0949e-13 | 2.7898e-12 | 1.2619e-11 |
| optimized_prepared | 9.0949e-13 | 2.7928e-12 | 1.2875e-11 |
| optimized_public | 9.0949e-13 | 2.7928e-12 | 1.2875e-11 |
| spline_approximation | 9.0949e-13 | 7.4482e-09 | 2.2911e-08 |

| Boundary | Median ms | Min–max ms | Bytes/call | Allocations | Samples |
|---|---:|---:|---:|---:|---:|
| ACE polynomial | 20.0264 | 17.8291–33.7733 | 64,313,216 | 183,925 | 10 |
| ETACE polynomial | 193.3885 | 181.6183–351.1577 | 205,989,144 | 124,145 | 10 |
| ACE spline | 12.8230 | 11.4258–13.7031 | 31,862,176 | 111,353 | 10 |
| PFW prefix public | 16.3077 | 16.2816–16.6651 | 4,952,632 | 92,736 | 10 |
| PFW DAG public | 13.3634 | 13.3175–13.8133 | 4,952,632 | 92,736 | 10 |
| PFW prefix prepared | 12.3415 | 12.2658–12.5356 | 32 | 1 | 10 |
| PFW DAG prepared | 9.3140 | 9.2893–9.6177 | 32 | 1 | 10 |

Setup elapsed seconds (excluded from warmed trials): `{'dag_calculator_s': 0.000496932, 'dag_executor_s': 0.000419452, 'optimized_calculator_s': 0.387453988, 'optimized_executor_s': 0.024145627}`.

## bcc-1024

Composition in Ti/Al/V order: [341, 340, 343]; directed neighbors: 58788. Face heights: [27.46290324036758, 27.0601365893845, 27.363200000000003] Å; minimum distance: 2.85925636738 Å. Verified center-1 V–V nearest image: atom 114, shift [0, -1, 0], distance 2.88172434261 Å.

Cell (column-major): `[27.471999999999998, 0.0, 0.0, 0.6255999999999999, 27.064, 0.0, -0.29919999999999997, 0.46240000000000003, 27.3632]`. Geometry SHA256: `d6bee510061d378412594ff00dce47e732a20a3dfff34a59692d311eb68f17b6`.

| Comparison | Energy/atom | Force component | Virial component |
|---|---:|---:|---:|
| ET_conversion | 1.3642e-12 | 1.3808e-14 | 1.0207e-11 |
| dag_prepared | 9.0949e-13 | 2.7414e-12 | 5.9948e-11 |
| dag_public | 9.0949e-13 | 2.7414e-12 | 5.9948e-11 |
| optimized_prepared | 9.0949e-13 | 2.736e-12 | 6.2144e-11 |
| optimized_public | 9.0949e-13 | 2.736e-12 | 6.2144e-11 |
| spline_approximation | 4.5475e-13 | 1.1383e-08 | 4.1279e-08 |

| Boundary | Median ms | Min–max ms | Bytes/call | Allocations | Samples |
|---|---:|---:|---:|---:|---:|
| ACE polynomial | 155.7826 | 154.1678–169.6844 | 516,592,656 | 1,563,115 | 10 |
| ETACE polynomial | 2166.3011 | 2076.8135–2239.6730 | 1,652,324,240 | 1,286,588 | 10 |
| ACE spline | 102.0470 | 91.5146–110.8728 | 256,470,224 | 981,379 | 10 |
| PFW prefix public | 134.1439 | 132.3844–136.8798 | 46,507,408 | 1,164,238 | 10 |
| PFW DAG public | 109.3799 | 108.2198–112.0964 | 46,507,408 | 1,164,238 | 10 |
| PFW prefix prepared | 98.3554 | 98.2115–98.5464 | 32 | 1 | 10 |
| PFW DAG prepared | 74.6292 | 74.4072–74.9667 | 32 | 1 | 10 |

Setup elapsed seconds (excluded from warmed trials): `{'dag_calculator_s': 0.000541192, 'dag_executor_s': 0.000481092, 'optimized_calculator_s': 0.000587052, 'optimized_executor_s': 0.000504622}`.

## bcc-3456

Composition in Ti/Al/V order: [1151, 1151, 1154]; directed neighbors: 197952. Face heights: [41.19435486055137, 40.59020488407674, 41.044799999999995] Å; minimum distance: 2.85645930622 Å. Verified center-1 V–V nearest image: atom 2, shift [0, 0, 0], distance 2.89633634636 Å.

Cell (column-major): `[41.208, 0.0, 0.0, 0.9383999999999999, 40.596, 0.0, -0.4487999999999999, 0.6936, 41.044799999999995]`. Geometry SHA256: `5426e0f29eb821da3d17069408b9d6714e8679d2339ff72977a74e2bfe5c9939`.

| Comparison | Energy/atom | Force component | Virial component |
|---|---:|---:|---:|
| ET_conversion | 2.6948e-13 | 1.6875e-14 | 8.7311e-11 |
| dag_prepared | 2.6948e-13 | 3.1941e-12 | 1.1278e-10 |
| dag_public | 2.6948e-13 | 3.1941e-12 | 1.1278e-10 |
| optimized_prepared | 2.6948e-13 | 3.1931e-12 | 1.1187e-10 |
| optimized_public | 2.6948e-13 | 3.1931e-12 | 1.1187e-10 |
| spline_approximation | 0 | 8.2202e-09 | 4.5315e-07 |

| Boundary | Median ms | Min–max ms | Bytes/call | Allocations | Samples |
|---|---:|---:|---:|---:|---:|
| ACE polynomial | 520.7956 | 518.2001–568.7844 | 1,743,957,008 | 5,478,182 | 10 |
| ETACE polynomial | 7394.7550 | 7078.9143–7431.8976 | 5,577,945,632 | 5,031,646 | 10 |
| ACE spline | 339.7365 | 327.8773–488.5606 | 868,071,440 | 3,519,398 | 10 |
| PFW prefix public | 461.5062 | 448.1533–476.5325 | 173,034,112 | 4,920,316 | 10 |
| PFW DAG public | 378.9566 | 365.9030–393.6133 | 173,034,112 | 4,920,316 | 10 |
| PFW prefix prepared | 328.0636 | 327.8665–329.1624 | 32 | 1 | 10 |
| PFW DAG prepared | 248.2504 | 247.9631–248.4179 | 32 | 1 | 10 |

Setup elapsed seconds (excluded from warmed trials): `{'dag_calculator_s': 0.000497842, 'dag_executor_s': 0.000499972, 'optimized_calculator_s': 0.000557712, 'optimized_executor_s': 0.000448552}`.

## hcp-150

Composition in Ti/Al/V order: [49, 49, 52]; directed neighbors: 5700. Face heights: [13.798675146155015, 13.785156303696333, 15.770794759935214] Å; minimum distance: 3.11909070676 Å. Verified center-1 V–V nearest image: atom 2, shift [0, 0, 0], distance 3.14824376125 Å.

Cell (column-major): `[16.16, 0.0, 0.0, 8.398697348592673, 13.787124428248262, 0.0, -0.17244407789193575, 0.2665044840148098, 15.770794759935216]`. Geometry SHA256: `83d4dd7121f8dd1f8597833981383a9da99df6d11e95e3ed04240ede6234af70`.

| Comparison | Energy/atom | Force component | Virial component |
|---|---:|---:|---:|
| ET_conversion | 1.9403e-13 | 8.6597e-15 | 2.2169e-12 |
| dag_prepared | 1.9403e-13 | 2.3324e-12 | 1.1312e-11 |
| dag_public | 1.9403e-13 | 2.3324e-12 | 1.1312e-11 |
| optimized_prepared | 1.9403e-13 | 2.3345e-12 | 1.1084e-11 |
| optimized_public | 1.9403e-13 | 2.3345e-12 | 1.1084e-11 |
| spline_approximation | 3.8805e-13 | 2.0932e-09 | 1.4462e-08 |

| Boundary | Median ms | Min–max ms | Bytes/call | Allocations | Samples |
|---|---:|---:|---:|---:|---:|
| ACE polynomial | 18.0110 | 15.9290–29.1417 | 56,039,424 | 148,917 | 10 |
| ETACE polynomial | 144.3116 | 139.0629–306.8069 | 171,643,144 | 99,035 | 10 |
| ACE spline | 10.8761 | 10.6560–12.6340 | 30,839,424 | 92,817 | 10 |
| PFW prefix public | 14.9337 | 14.9053–15.2444 | 4,322,456 | 76,901 | 10 |
| PFW DAG public | 11.4064 | 11.3736–11.6255 | 4,322,456 | 76,901 | 10 |
| PFW prefix prepared | 11.6430 | 11.6150–11.8411 | 32 | 1 | 10 |
| PFW DAG prepared | 8.1182 | 8.0879–8.2742 | 32 | 1 | 10 |

Setup elapsed seconds (excluded from warmed trials): `{'dag_calculator_s': 0.000435232, 'dag_executor_s': 0.000421202, 'optimized_calculator_s': 0.000514852, 'optimized_executor_s': 0.000434682}`.

## omega-240

Composition in Ti/Al/V order: [78, 79, 83]; directed neighbors: 13272. Face heights: [15.868476418078265, 15.852929749250785, 14.3355] Å; minimum distance: 2.60346346371 Å. Verified center-1 V–V nearest image: atom 49, shift [0, 0, 0], distance 2.8213065704 Å.

Cell (column-major): `[18.584, 0.0, 0.0, 9.658501950881574, 15.855193092485502, 0.0, -0.15675, 0.24225000000000002, 14.3355]`. Geometry SHA256: `d723432754d6cbfca2d031dc0d9f3689ca996d55866f2cf7e89707cb63d27611`.

| Comparison | Energy/atom | Force component | Virial component |
|---|---:|---:|---:|
| ET_conversion | 2.4253e-13 | 1.5168e-14 | 1.3642e-12 |
| dag_prepared | 2.4253e-13 | 3.0565e-12 | 1.5916e-11 |
| dag_public | 2.4253e-13 | 3.0565e-12 | 1.5916e-11 |
| optimized_prepared | 2.4253e-13 | 3.0575e-12 | 1.603e-11 |
| optimized_public | 2.4253e-13 | 3.0575e-12 | 1.603e-11 |
| spline_approximation | 7.276e-13 | 4.345e-09 | 4.645e-08 |

| Boundary | Median ms | Min–max ms | Bytes/call | Allocations | Samples |
|---|---:|---:|---:|---:|---:|
| ACE polynomial | 38.8086 | 35.1718–46.5150 | 117,802,944 | 333,603 | 10 |
| ETACE polynomial | 350.4543 | 330.4416–522.4800 | 376,972,504 | 224,417 | 10 |
| ACE spline | 21.7375 | 21.1666–29.2258 | 59,081,280 | 202,323 | 10 |
| PFW prefix public | 30.2214 | 30.1097–30.4546 | 9,076,776 | 168,323 | 10 |
| PFW DAG public | 24.4937 | 24.3816–24.8133 | 9,076,776 | 168,323 | 10 |
| PFW prefix prepared | 22.7644 | 22.7389–22.8961 | 32 | 1 | 10 |
| PFW DAG prepared | 17.1317 | 17.0444–17.3817 | 32 | 1 | 10 |

Setup elapsed seconds (excluded from warmed trials): `{'dag_calculator_s': 0.000470482, 'dag_executor_s': 0.000454912, 'optimized_calculator_s': 0.000545222, 'optimized_executor_s': 0.000440882}`.

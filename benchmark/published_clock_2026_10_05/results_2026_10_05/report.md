# Published-clock benchmark: shared-record comparison

Direction agreement: published-significant contrasts (T1) or expected direction (T2). Call balanced accuracy: mean recall over the published calls present (higher, lower, or none = not significant), each system calling its direction when significant and none otherwise; S0 (constant higher, significant) scores at most 1/k with k classes present. Significance agreement applies the paper's stated threshold, else P < 0.05 (count of assumed thresholds shown). Brackets: 95% bootstrap intervals over papers (T1) or series (T2); (n) = denominator.

## T1, all records

| System | Covered / records | Direction agreement | Call balanced accuracy (higher/lower/none n) | Significance agreement | Assumed thresholds |
| --- | --- | --- | --- | --- | --- |
| S0 | 580 / 580 | 0.764 [0.644, 0.873] (225) | 0.333 [0.333, 0.333] (172/53/260) | 0.464 (485) | 0 |
| S1 | 107 / 580 | 0.857 [0.704, 0.979] (35) | 0.669 [0.526, 0.828] (25/10/54) | 0.809 (89) | 71 |
| S2 | 116 / 580 | 0.875 [0.744, 0.965] (40) | 0.578 [0.431, 0.733] (30/10/54) | 0.702 (94) | 74 |
| S3 | 154 / 580 | 0.963 [0.914, 1.000] (54) | 0.673 [0.530, 0.798] (38/16/59) | 0.770 (113) | 89 |
| S4:openai:gpt-4.1 | 323 / 580 | 0.767 [0.615, 0.883] (103) | 0.553 [0.473, 0.649] (85/18/154) | 0.634 (257) | 0 |
| S4:openai:gpt-5.5 | 324 / 580 | 0.683 [0.500, 0.833] (104) | 0.489 [0.371, 0.610] (85/19/154) | 0.620 (258) | 0 |
| S4:openai:gpt-6-astra | 324 / 580 | 0.673 [0.525, 0.793] (104) | 0.532 [0.417, 0.650] (85/19/154) | 0.694 (258) | 0 |
| S4:openai:gpt-6-luna | 324 / 580 | 0.567 [0.382, 0.725] (104) | 0.461 [0.358, 0.572] (85/19/154) | 0.632 (258) | 0 |
| S4:openai:gpt-6-sol | 324 / 580 | 0.673 [0.500, 0.813] (104) | 0.494 [0.383, 0.611] (85/19/154) | 0.663 (258) | 0 |

## T1, records with a ClockBase clock

| System | Covered / records | Direction agreement | Call balanced accuracy (higher/lower/none n) | Significance agreement | Assumed thresholds |
| --- | --- | --- | --- | --- | --- |
| S0 | 327 / 327 | 0.815 [0.676, 0.915] (108) | 0.333 [0.333, 0.333] (88/20/158) | 0.406 (266) | 0 |
| S1 | 107 / 327 | 0.857 [0.703, 0.976] (35) | 0.669 [0.527, 0.819] (25/10/54) | 0.809 (89) | 71 |
| S2 | 116 / 327 | 0.875 [0.750, 0.971] (40) | 0.578 [0.439, 0.731] (30/10/54) | 0.702 (94) | 74 |
| S3 | 154 / 327 | 0.963 [0.913, 1.000] (54) | 0.673 [0.533, 0.799] (38/16/59) | 0.770 (113) | 89 |
| S4:openai:gpt-4.1 | 301 / 327 | 0.745 [0.587, 0.864] (94) | 0.540 [0.454, 0.639] (76/18/149) | 0.617 (243) | 0 |
| S4:openai:gpt-5.5 | 302 / 327 | 0.674 [0.477, 0.816] (95) | 0.485 [0.367, 0.608] (76/19/149) | 0.623 (244) | 0 |
| S4:openai:gpt-6-astra | 302 / 327 | 0.653 [0.483, 0.782] (95) | 0.521 [0.407, 0.643] (76/19/149) | 0.689 (244) | 0 |
| S4:openai:gpt-6-luna | 302 / 327 | 0.547 [0.349, 0.707] (95) | 0.459 [0.349, 0.569] (76/19/149) | 0.639 (244) | 0 |
| S4:openai:gpt-6-sol | 302 / 327 | 0.653 [0.454, 0.783] (95) | 0.492 [0.374, 0.615] (76/19/149) | 0.668 (244) | 0 |

## T1, records covered by every system of S1, S2, S3 (60 records, 26 papers)

| System | Direction agreement | Call balanced accuracy (higher/lower/none n) | Significance agreement |
| --- | --- | --- | --- |
| S0 | 0.680 [0.409, 0.889] (25) | 0.333 [0.333, 0.333] (17/8/24) | 0.510 (49) |
| S1 | 0.880 [0.700, 1.000] (25) | 0.708 [0.527, 0.867] (17/8/24) | 0.776 (49) |
| S2 | 0.960 [0.863, 1.000] (25) | 0.632 [0.433, 0.806] (17/8/24) | 0.673 (49) |
| S3 | 1.000 [1.000, 1.000] (25) | 0.694 [0.511, 0.854] (17/8/24) | 0.755 (49) |

## T1, records covered by every system of S3, S4:openai:gpt-4.1, S4:openai:gpt-5.5, S4:openai:gpt-6-astra, S4:openai:gpt-6-luna, S4:openai:gpt-6-sol (154 records, 41 papers)

| System | Direction agreement | Call balanced accuracy (higher/lower/none n) | Significance agreement |
| --- | --- | --- | --- |
| S0 | 0.704 [0.488, 0.857] (54) | 0.333 [0.333, 0.333] (38/16/61) | 0.470 (115) |
| S3 | 0.963 [0.917, 1.000] (54) | 0.673 [0.527, 0.805] (38/16/59) | 0.770 (113) |
| S4:openai:gpt-4.1 | 0.704 [0.522, 0.848] (54) | 0.516 [0.416, 0.616] (38/16/61) | 0.557 (115) |
| S4:openai:gpt-5.5 | 0.574 [0.333, 0.771] (54) | 0.507 [0.366, 0.645] (38/16/61) | 0.600 (115) |
| S4:openai:gpt-6-astra | 0.556 [0.370, 0.700] (54) | 0.518 [0.384, 0.651] (38/16/61) | 0.661 (115) |
| S4:openai:gpt-6-luna | 0.481 [0.281, 0.648] (54) | 0.455 [0.327, 0.574] (38/16/61) | 0.583 (115) |
| S4:openai:gpt-6-sol | 0.611 [0.397, 0.769] (54) | 0.496 [0.356, 0.622] (38/16/61) | 0.626 (115) |

## T1, records covered by every system of S1, S2, S3, S4:openai:gpt-4.1, S4:openai:gpt-5.5, S4:openai:gpt-6-astra, S4:openai:gpt-6-luna, S4:openai:gpt-6-sol (60 records, 26 papers)

| System | Direction agreement | Call balanced accuracy (higher/lower/none n) | Significance agreement |
| --- | --- | --- | --- |
| S0 | 0.680 [0.409, 0.889] (25) | 0.333 [0.333, 0.333] (17/8/24) | 0.510 (49) |
| S1 | 0.880 [0.700, 1.000] (25) | 0.708 [0.527, 0.867] (17/8/24) | 0.776 (49) |
| S2 | 0.960 [0.863, 1.000] (25) | 0.632 [0.433, 0.806] (17/8/24) | 0.673 (49) |
| S3 | 1.000 [1.000, 1.000] (25) | 0.694 [0.511, 0.854] (17/8/24) | 0.755 (49) |
| S4:openai:gpt-4.1 | 0.680 [0.412, 0.900] (25) | 0.511 [0.382, 0.651] (17/8/24) | 0.612 (49) |
| S4:openai:gpt-5.5 | 0.640 [0.381, 0.849] (25) | 0.551 [0.387, 0.719] (17/8/24) | 0.673 (49) |
| S4:openai:gpt-6-astra | 0.600 [0.333, 0.810] (25) | 0.587 [0.417, 0.751] (17/8/24) | 0.735 (49) |
| S4:openai:gpt-6-luna | 0.560 [0.316, 0.769] (25) | 0.512 [0.352, 0.645] (17/8/24) | 0.633 (49) |
| S4:openai:gpt-6-sol | 0.640 [0.368, 0.852] (25) | 0.599 [0.433, 0.764] (17/8/24) | 0.735 (49) |

## T1, LLM-only predictions by whether the model said it recognized the publication

| Model and recognition | Covered | Direction agreement | Call balanced accuracy | Significance agreement |
| --- | --- | --- | --- | --- |
| S4:openai:gpt-4.1|recognized=True | 50 | 0.958 [0.750, 1.000] (24) | 0.600 [0.333, 0.667] (19/5/6) | 0.800 (30) |
| S4:openai:gpt-4.1|recognized=False | 273 | 0.709 [0.526, 0.842] (79) | 0.477 [0.405, 0.547] (66/13/148) | 0.612 (227) |
| S4:openai:gpt-5.5|recognized=True | 66 | 1.000 [1.000, 1.000] (29) | 0.720 [0.485, 0.823] (24/5/19) | 0.646 (48) |
| S4:openai:gpt-5.5|recognized=False | 258 | 0.560 [0.329, 0.747] (75) | 0.375 [0.295, 0.482] (61/14/135) | 0.614 (210) |
| S4:openai:gpt-6-astra|recognized=True | 65 | 0.933 [0.891, 1.000] (30) | 0.849 [0.737, 0.906] (25/5/17) | 0.830 (47) |
| S4:openai:gpt-6-astra|recognized=False | 259 | 0.568 [0.358, 0.733] (74) | 0.400 [0.332, 0.485] (60/14/137) | 0.664 (211) |
| S4:openai:gpt-6-luna|recognized=True | 42 | 0.812 [0.357, 1.000] (16) | 0.667 [0.417, 0.867] (12/4/8) | 0.625 (24) |
| S4:openai:gpt-6-luna|recognized=False | 282 | 0.523 [0.324, 0.696] (88) | 0.381 [0.319, 0.456] (73/15/146) | 0.632 (234) |
| S4:openai:gpt-6-sol|recognized=True | 50 | 0.962 [0.909, 1.000] (26) | 0.784 [0.667, 0.889] (21/5/10) | 0.833 (36) |
| S4:openai:gpt-6-sol|recognized=False | 274 | 0.577 [0.333, 0.764] (78) | 0.355 [0.303, 0.419] (64/14/144) | 0.635 (222) |

## T1, S3 against each model on records the model did not recognize (and S3 covers)

| Model | Records (papers) | System | Direction agreement | Call balanced accuracy | Significance agreement |
| --- | --- | --- | --- | --- | --- |
| S4:openai:gpt-4.1 | 123 (36) | S0 | 0.725 [0.481, 0.902] (40) | 0.333 [0.333, 0.333] (29/11/57) | 0.412 (97) |
| S4:openai:gpt-4.1 | 123 (36) | S3 | 0.975 [0.900, 1.000] (40) | 0.697 [0.529, 0.835] (29/11/55) | 0.779 (95) |
| S4:openai:gpt-4.1 | 123 (36) | S4:openai:gpt-4.1 | 0.625 [0.421, 0.820] (40) | 0.437 [0.356, 0.526] (29/11/57) | 0.515 (97) |
| S4:openai:gpt-5.5 | 116 (34) | S0 | 0.711 [0.438, 0.894] (38) | 0.333 [0.333, 0.333] (27/11/52) | 0.422 (90) |
| S4:openai:gpt-5.5 | 116 (34) | S3 | 0.974 [0.895, 1.000] (38) | 0.674 [0.497, 0.821] (27/11/50) | 0.773 (88) |
| S4:openai:gpt-5.5 | 116 (34) | S4:openai:gpt-5.5 | 0.395 [0.138, 0.625] (38) | 0.379 [0.265, 0.498] (27/11/52) | 0.578 (90) |
| S4:openai:gpt-6-astra | 117 (35) | S0 | 0.711 [0.440, 0.896] (38) | 0.333 [0.333, 0.333] (27/11/53) | 0.418 (91) |
| S4:openai:gpt-6-astra | 117 (35) | S3 | 0.974 [0.905, 1.000] (38) | 0.668 [0.491, 0.810] (27/11/51) | 0.764 (89) |
| S4:openai:gpt-6-astra | 117 (35) | S4:openai:gpt-6-astra | 0.395 [0.190, 0.559] (38) | 0.369 [0.296, 0.433] (27/11/53) | 0.615 (91) |
| S4:openai:gpt-6-luna | 124 (37) | S0 | 0.707 [0.500, 0.875] (41) | 0.333 [0.333, 0.333] (29/12/56) | 0.423 (97) |
| S4:openai:gpt-6-luna | 124 (37) | S3 | 0.951 [0.870, 1.000] (41) | 0.681 [0.519, 0.824] (29/12/54) | 0.779 (95) |
| S4:openai:gpt-6-luna | 124 (37) | S4:openai:gpt-6-luna | 0.366 [0.179, 0.531] (41) | 0.353 [0.266, 0.443] (29/12/56) | 0.577 (97) |
| S4:openai:gpt-6-sol | 124 (36) | S0 | 0.718 [0.458, 0.903] (39) | 0.333 [0.333, 0.333] (28/11/56) | 0.411 (95) |
| S4:openai:gpt-6-sol | 124 (36) | S3 | 0.974 [0.893, 1.000] (39) | 0.674 [0.499, 0.817] (28/11/54) | 0.774 (93) |
| S4:openai:gpt-6-sol | 124 (36) | S4:openai:gpt-6-sol | 0.487 [0.214, 0.694] (39) | 0.351 [0.258, 0.436] (28/11/56) | 0.589 (95) |

## T2 (ComputAgeBench expected directions, one record per ClockBase clock)

| System | Covered / records | Expected direction | Significant in expected direction |
| --- | --- | --- | --- |
| S1 | 73 / 157 | 0.808 [0.661, 0.929] (73) | 0.370 [0.210, 0.562] |
| S2 | 72 / 157 | 0.722 [0.569, 0.875] (72) | 0.347 [0.200, 0.525] |
| S3 | 70 / 157 | 0.800 [0.622, 0.924] (70) | 0.343 [0.212, 0.522] |


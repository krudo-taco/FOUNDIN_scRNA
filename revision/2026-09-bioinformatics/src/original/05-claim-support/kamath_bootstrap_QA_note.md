# Verification of two secondary bootstrap-floor results

This diagnostic was added after viewing the frozen Kamath analysis. The prespecified results remain unchanged in `results/kamath_all_edges.csv`. It does not create a replacement primary analysis.

The Sepulveda-only DA-minus-NonDA comparisons for TOLLIP (RNA/QC/library adjustment) and ATG4B (rank/QC/library adjustment) have zero tail exceedances in 19,999 Rademacher wild-bootstrap draws. The reported add-one P value is therefore 1/20,000 for each, with BH q = 0.0252 across the 1,008 planned edges. Independent computation reproduces these values and finds that all 2,048 possible sign patterns were already represented in the Monte Carlo sample.

Exhaustive evaluation of all 2^11 sign patterns also produces zero exceedances. TOLLIP's observed absolute t is 4.472385, just above the largest reference absolute t of 4.459440; ATG4B's corresponding values are 2.303181 and 2.273200. Thus the extreme Monte Carlo P values do not result from missing a rare sign pattern, and are not a transcription or programming error. They lie at an add-one floor of a finite reference distribution.

Using the add-one convention with all 2,048 unique patterns gives 1/2,049 = 0.000488 for each. In a bounded floor sensitivity that substitutes only these two values and leaves other P values unchanged, both BH values become 0.168 and no edge reaches BH < 0.05. This calculation demonstrates sensitivity to the floor convention; it is not an exact disease-label randomization test or a new, fully enumerated primary analysis.

The HC3 t-based intervals and P values provide additional context. TOLLIP has coefficient 0.060871, 95% interval 0.028688 to 0.093055, and nominal P = 0.002893. ATG4B has coefficient 0.080102, interval -0.002137 to 0.162341, and nominal P = 0.054734. Neither gene is in the original fixed 16-gene direction profile. These are secondary, method-sensitive observations and cannot establish independent replication of the original PD-specific, DA-specific coupling claim.

All four primary RNA panel/profile tests, and all 48 planned panel/profile tests, remain nonsignificant after their specified corrections. `results/kamath_all_edges_with_QA_flags.csv` retains the original estimates, P values and corrections alongside the diagnostic flags and bounded floor sensitivity. No cell, donor, gene or model was changed to obtain a desired outcome.

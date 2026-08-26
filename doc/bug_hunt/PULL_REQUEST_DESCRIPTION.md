# Triage and fixes for the oldest open issues

This branch processes the **150 oldest open issues** of the repository
(three consecutive batches of 50, from #998 up to #2537, spanning
January 2019 -- May 2024). Each issue body was read from the GitHub API,
ambiguous cases were checked against their comment threads, and every
classification was cross-checked against the current source tree.

Each issue received:

- a **category**: actual *bug*, *feature request*, or *other*
  (documentation work, usage questions, duplicates);
- an **effort estimate**: `effort:easy`, `effort:medium`,
  `effort:hard`, or *no longer valid*;
- where useful, an **area** qualifier (`probas`, `surrogate`,
  `sensitivity`, `graph`, `API`, `process`, `example`, `theory`, ...).

Thirty contained issues were fixed on this branch, one commit each,
with Python test cases added for every behavioural change and the
affected `.expout` references regenerated when required. The full
classification details are available in
[triage_oldest_50_issues.pdf](./triage_oldest_50_issues.pdf).

---

## 1. Issues fixed (one commit per issue)

### Bugs

| Issue | Fix |
|---|---|
| #1357 | QQ/PP/CDF plots used an almost invisible marker for large samples; outliers are now visible (`plus` style above `VisualTest-CloudMediumSize`) |
| #1479 | Nested `CompositeDistribution` composes both functions so `getAntecedent()` always returns a non-composite law |
| #1596 | Mixed mixtures (continuous + discrete) can be drawn: continuous-part curve plus one arrow per support point of every discrete component, normalized scale, infinite supports handled |
| #1678 | The process Monte-Carlo example now builds the bidimensionnal white noise described by its introduction, so the exact probability 0.17008 is actually comparable to the estimate |
| #1841 | AbdoRackwitz stopping rule documented (disjunction of two criteria pairs); convergence status message names the pair that fired; iteration exhaustion reported explicitly |
| #2432 | Adding/subtracting `Point`s of different dimensions keeps the informative dimension-mismatch message instead of a generic "unsupported operand" TypeError |
| #2498 | Sample-versus-model QQ plots no longer drop the largest observation (probability capped at 1 - machine epsilon) |
| #2502 | Out-of-bounds access `(i, j)` on symmetric matrices reports the offending column index instead of the swapped internal row index |

### Enhancements

| Issue | Fix |
|---|---|
| #1243 | `UserDefinedStationaryCovarianceModel(mesh, covariance)` constructor taking a single matrix replicated over the mesh |
| #1363 | Optional `yLabel` argument for `DrawCorrelationCoefficients` so SRC indices are not plotted as "correlation coefficient" |
| #1424 | `ProcessSample::split(index)`, mirroring `Sample::split` |
| #1463 | `Field::getDimension()` added; `getInputDimension`/`getOutputDimension` deprecated |
| #1475 | `Field +/- Field` returns a proper field on the common mesh, with mesh and dimension validation |
| #1489 | `Mixture::drawPDF` override handling components of both natures (curve + vertical arrows head-up for atoms) |
| #1554 | `KernelSmoothing.buildWeighted(sample, weights)` for non parametric importance sampling; bandwidth selection stays unweighted by design |
| #1583 | `Graph.setXMin/setXMax/setYMin/setYMax` convenience bounds accessors |
| #1730 | `Sample.erase(indices)` / `ProcessSample.erase(indices)` removing a set of positions in one call |
| #1734 | `computeMinimumVolumeLevelSetCollectionWithThreshold(probs)` computing several minimum-volume level sets at once, sharing one QMC sample between levels |
| #1937 | `FunctionalChaosSobolIndices.draw()` drawing first/total order indices through `DrawSobolIndices` |
| #1996 | `FunctionCollection` slicing and sequence indexing (`collection[1:3]`, negative and stepped slices, index lists) |
| #2119 | `ComposedFunction.get/setLeftFunction` and `get/setRightFunction` |

### Documentation

| Issue | Fix |
|---|---|
| #998 | New theory page on piecewise linear, cubic Hermite and P1 Lagrange interpolation |
| #1318 | Gallery example: functions with discrete inputs and the finite-difference pitfall |
| #1322 | Notes on initial-versus-transformed data estimation (moments vs least squares) in docstrings and theory page |
| #1410 | SphericalModel help page documents positive definiteness limits (Stein 1999) |
| #1477 | KissFFT usage example with doctest-verified outputs |
| #1496 | EfficientGlobalOptimization example showing multistart tuning and solver replacement |
| #1511 | LeastSquaresMethod docstrings state the weighted problem actually solved and use a consistent notation |
| #1582 | Stochastic process definitions: spatial mean renamed, inconsistent `g` notation removed, reference added |
| #1707 | Warning added: ANCOVA indices are not restricted to [0, 1] nor summing to one with correlated inputs |

A follow-up commit also restores the `drawPDF` overloads that SWIG
shadowed on `Mixture` once the new override was declared.

---

## 2. Invalid or already resolved issues (proposed closure)

| Issue | Reason |
|---|---|
| #1077 | Constructor docstring notation already rewritten upstream (`(*spatialDim*)` form) |
| #1231 | Report predates major renamings; uses pre-1.14 API names, not reproducible |
| #1358 | `CorrelationAnalysis` has been documented upstream in the meantime |
| #1396 | Duplicate of another report |
| #1427 | Index-based extraction of `Sample`/`ProcessSample` already works in master |
| #1245 | Comparison operator fix landed in another pull request |
| #1468 | Usage question answered in the thread |
| #1542 | Support question; addressed by later HMAT / Karhunen-Loeve developments |
| #1585 | Superseded: p-values now estimated by the adaptive Lilliefors scheme (`FittingTest-Lilliefors*` keys) |
| #2521 | Verified fixed in master: `Normal::computeProbability` handles infinite bounds in 1D |
| #2537 | Verified fixed in master: discrete distributions accepted through the factory signature of `FittingTest.ChiSquared` |

---

## 3. Other issues examined but not fixed - proposed requalification

Proposed qualifiers follow the `type . area . effort:*` scheme; items
marked *(parked)* touch the Gaussian process / kriging area which is
evolving in parallel pull requests, and items marked *(policy)* require
an upstream API decision before any change.

### Batch 1 (#998 - #1585)

| Issue | Proposed qualifiers |
|---|---|
| #1056 | enhancement, sensitivity, effort:hard |
| #1225 | enhancement, surrogate, effort:hard (parked) |
| #1227 | enhancement, probas, effort:hard |
| #1230 | enhancement, API, python, effort:medium |
| #1237 | enhancement, probas, effort:hard |
| #1239 | enhancement, API, probas, effort:medium (design decision needed) |
| #1242 | bug, probas, refactor, effort:hard |
| #1244 | bug, performance, probas, effort:hard |
| #1246 | enhancement, probas, API, effort:hard |
| #1250 | bug, reliability, API, effort:hard |
| #1256 | enhancement, API, refactor, effort:hard |
| #1262 | enhancement, performance, effort:hard |
| #1284 | bug, probas, effort:medium (estimator claimed upstream) |
| #1320 | enhancement, surrogate, python, effort:hard |
| #1325 | enhancement, probas, API, effort:hard |
| #1328 | question, surrogate, needs review, effort:medium |
| #1353 | enhancement, surrogate, effort:medium (parked) |
| #1409 | enhancement, API, surrogate, effort:medium (API break) |
| #1445 | enhancement, probas, effort:hard |
| #1480 | doc, probas, API, effort:medium |
| #1481 | doc, optimization, effort:medium |
| #1488 | doc, example, effort:hard |

### Batch 2 (#1587 - #2164)

| Issue | Proposed qualifiers |
|---|---|
| #1587 | enhancement, graph, effort:medium |
| #1595 | doc, API, effort:hard |
| #1601 | enhancement, probas, effort:hard (changes all QMC outputs) |
| #1700 | question, surrogate, effort:medium (parked) |
| #1703 | bug, performance, process, effort:hard (parked) |
| #1720 | bug, probas, effort:hard |
| #1722 | doc, example, effort:medium |
| #1741 | doc, example, effort:medium |
| #1757 | doc, example, effort:medium |
| #1762 | doc, sensitivity, effort:easy |
| #1765 | enhancement, performance, effort:hard |
| #1769 | enhancement, graph, API, effort:medium |
| #1779 | enhancement, graph, effort:medium |
| #1780 | doc, example, python, effort:medium |
| #1786 | enhancement, surrogate, effort:medium (parked) |
| #1850 | enhancement, probas, effort:hard |
| #1852 | enhancement, surrogate, effort:medium |
| #1873 | bug, surrogate, performance, effort:hard (parked) |
| #1875 | enhancement, python, probas, effort:medium |
| #1886 | enhancement, graph, effort:medium |
| #1888 | enhancement, API, python, effort:hard |
| #1892 | doc, example, sensitivity, effort:easy |
| #1893 | enhancement, reliability, effort:medium |
| #1899 | doc, API, effort:medium |
| #1924 | doc, question, effort:medium |
| #1928 | doc, effort:easy |
| #1951 | doc, surrogate, effort:easy |
| #1957 | refactor, API, effort:medium |
| #1960 | doc, example, effort:medium |
| #1992 | doc, theory, effort:hard |
| #2000 | doc, effort:easy |
| #2078 | bug, probas, effort:hard |
| #2080 | doc, surrogate, effort:easy |
| #2114 | doc, effort:hard |
| #2118 | enhancement, graph, effort:medium |
| #2127 | doc, effort:medium (site configuration) |
| #2128 | refactor, API, effort:medium (renames need deprecation cycle) |
| #2164 | refactor, API, effort:medium (deprecation policy) |

### Batch 3 (#2178 - #2537)

| Issue | Proposed qualifiers |
|---|---|
| #2178 | doc, example, signal, effort:medium |
| #2201 | enhancement, API, python, effort:easy |
| #2228 | refactor, API, effort:medium (policy) |
| #2247 | doc, reliability, effort:easy |
| #2248 | doc, API, probas, effort:easy |
| #2249 | enhancement, optimization, effort:hard |
| #2253 | doc, example, surrogate, effort:easy |
| #2260 | doc, sensitivity, effort:medium |
| #2280 | doc, calibration, graph, effort:medium |
| #2294 | doc, theory, surrogate, effort:hard |
| #2319 | bug, surrogate, effort:hard (parked) |
| #2331 | doc, sensitivity, effort:easy |
| #2333 | enhancement, sensitivity, surrogate, effort:medium |
| #2334 | tests, surrogate, effort:easy |
| #2335 | doc, surrogate, effort:easy |
| #2355 | doc, API, surrogate, effort:medium |
| #2365 | doc, surrogate, effort:medium (parked) |
| #2370 | doc, effort:medium (website front page) |
| #2376 | enhancement, effort:hard |
| #2378 | enhancement, effort:medium |
| #2379 | doc, sensitivity, effort:easy |
| #2380 | enhancement, sensitivity, surrogate, effort:medium |
| #2383 | doc, example, effort:medium |
| #2385 | doc, surrogate, effort:medium |
| #2390 | enhancement, probas, effort:medium |
| #2400 | question, surrogate, effort:medium |
| #2401 | bug, probas, effort:medium |
| #2407 | enhancement, python, effort:medium |
| #2408 | doc, surrogate, effort:easy |
| #2409 | doc, surrogate, effort:easy |
| #2428 | bug, graph, effort:medium |
| #2436 | enhancement, effort:hard (HERE macro / traceback policy) |
| #2440 | doc, probas, effort:easy |
| #2444 | enhancement, probas, effort:hard |
| #2445 | bug, graph, difficulty:medium |
| #2465 | cleanup, doc, effort:easy (remove versus document) |
| #2482 | enhancement, graph, API, effort:medium (policy) |
| #2483 | doc, example, graph, effort:easy |
| #2484 | refactor, API, effort:medium (policy) |
| #2485 | enhancement, sample, effort:easy |
| #2494 | doc, effort:easy |
| #2509 | enhancement, probas, API, effort:easy |
| #2512 | doc, build, effort:easy (documentation infrastructure) |
| #2533 | refactor, reliability, API, effort:hard |

---

## Testing

- Full build (library + SWIG regeneration) after every issue.
- New or extended Python test cases for each behavioural fix
  (`ott.assert_almost_equal` assertions; `.expout` files regenerated
  where the printed output legitimately changed, e.g. QQ plots).
- Corresponding C++ test suites kept passing
  (`t_VisualTest_std`, `t_ComposedDistribution_std`, ...).
- Documentation changes validated by the Sphinx build (warnings treated
  as errors) including gallery execution and doctests.
- Python sources lint-free (`flake8`).

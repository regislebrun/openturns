LinearCombinationDistribution analysis report
=============================================

Branch ``fix_linear_combination_distribution``
Report generated on 2026-08-25.

Context
-------

The numerical fixes of the class (commit "LinearCombinationDistribution:
fix 5 bugs") and the conversion of the Python tests to assert based checks
surfaced several questions about the behaviour of the class. This report
explores five of them empirically. Every section below is backed by a
script stored next to this file; all scripts were run with the branch
build (Release, default resource map unless stated otherwise).

.. contents::


1. Discrete fusion: exact equality merging of supports
------------------------------------------------------

**Question.** The fusion of discrete atoms convolves their supports scaled
by their weights, then merges equal support points with an exact ``==``
comparison. Can mathematically identical points produced through different
accumulation paths differ by one ulp and remain duplicated?

**Method.** :file:`p1_fusion.py` fuses families of Binomial/Poisson atoms
with integer, rational and random fractional weights (40 randomized cases)
and compares the stored support size against the cardinality of the exact
support (enumeration rounded at 9 digits).

**Results.**

- No case of inflation was found: the stored support always matched the
  exact cardinality. Scaled supports are of the form :math:`k \times w`
  for integer :math:`k`, so duplicate paths are built from single
  multiplications which agree bitwise; the subsequent sequential sums then
  also agree.
- A different effect was uncovered instead: the *pre-merges* performed by
  ``poissonMap``/``binomialMap`` (same weight, or same couple (weight, p))
  combine the distributions *before* their supports get truncated. For
  example ``Poi(2) + Poi(2) + 2 * Poi(2)`` first becomes
  ``Poi(4) + 2 * Poi(2)``, whose fused support holds 65 points (maximum
  64), while the convolution of the three originally truncated supports
  would reach 76 with 77 points. The lost tail mass is at the level of the
  per-atom truncation (~:math:`10^{-16}` here), so probabilities stay
  consistent (sum = 0.99999999999999991).

**Conclusion / recommendation.** No change needed for the ``==`` merging.
The pre-merge/truncation interplay is benign (it can only shrink the
fused support, which in turn favours further fusion) but its semantics
should be kept in mind if ``MaximumSupportSize`` is tuned tightly.

*Post-scriptum.* The remaining gap between the 180 raw convolution sums and
the 130 stored points comes from
:meth:`FiniteDiscreteDistribution.compactSupport`, called automatically by
the constructor in dimension 1: near-duplicate points within
``Distribution-SupportEpsilon`` are merged with their probabilities.
A tolerance based merging already exists one level below the fusion.


2. Range truncation versus heavy-tailed atoms
---------------------------------------------

**Question.** The range is clipped to
[position indicator :math:`\pm` beta * dispersion indicator] even when the
exact range is unbounded: beyond it the CDF returns exactly 1 and the PDF
exactly 0. How much mass is silently discarded for heavy-tailed atoms?

**Method.** :file:`p2_tails.py` evaluates the discarded mass on 500000
samples and checks whether :meth:`computeQuantile`\ ``(1 - 1e-3)`` hits the
clipped bound.

**Results.**

==========================  ========================  =======================  =====================
atoms                       upper bound               MC mass above            quantile(0.999)
==========================  ========================  =======================  =====================
LogNormal(sigmaLog=1.0)+N   20.04                     1.37e-3                  20.0386 (at bound)
LogNormal(sigmaLog=1.5)+N   79.36                     1.77e-3                  79.3624 (at bound)
LogNormal(sigmaLog=2.0)+N   467.2                     1.03e-3                  467.204 (at bound)
LogNormal(sigmaLog=2.5)+N   4422                      3.7e-4                   2265.57
Pareto(alpha=2.5)+N         14.37                     1.25e-3                  14.3662 (at bound)
==========================  ========================  =======================  =====================

**Conclusion / recommendation.** For mixtures with heavy-tailed atoms,
up to ~2e-3 of probability mass is declared zero beyond the clipped range,
and extreme quantiles (:math:`p \gtrsim 0.998`) collapse onto the bound.
This is 7 to 8 orders of magnitude above ``DefaultCDFEpsilon``.
Options: enlarge beta for subexponential atom families, expose the discarded
mass, or document that extreme quantiles of such mixtures require
``setBeta``/``setAlpha`` tuning.


3. Characteristic function cache economics
------------------------------------------

**Question.** The cache of characteristic function values stores whole
levels of the grid; its cumulative size after level m is about
:math:`(2m+1)^d`. With the default ``DefaultMaxSize=65536`` the fully
cached horizon is m <= 127 in dimension 2 but only **m <= 19 in dimension 3**.
What does this cost?

**Method.** :file:`p3_cache.py` times pointwise PDF evaluations of the
reference 4-atom 3D mixture (the grid3d test case), whose evaluation stops
on the ``MaximumPDFLevel`` bound after tens of levels, with the default
cache and with a 16M cache.

**Results.**

===============  =====================  ======================
maxSize          5 distinct points      20 distinct points
===============  =====================  ======================
65536 (default)  18.25 s (3651 ms/call) 73.18 s (3659 ms/call)
16777216         2.58 s (516 ms/call)   3.52 s (176 ms/call)
===============  =====================  ======================

Beyond the horizon every call recomputes levels ~20 to ~64 (about
:math:`10^6` characteristic evaluations): the observed speed-up of the
larger cache ranges from 7x to 20x, and warm calls get cheaper as the deep
levels get cached once.

**Conclusion / recommendation.** In dimension 3 the default cache is far
too small to cover even the bounded pointwise evaluation. Options: scale
``DefaultMaxSize`` with the dimension (e.g. document
:math:`(2 b_{max})^{dim}` as a rule of thumb), let the cache grow up to a
memory budget, or warn when the horizon is exceeded during an evaluation.


4. The SmallSize threshold applies to computeProbability only
-------------------------------------------------------------

**Question.** For 1D collections larger than ``SmallSize`` (100),
:meth:`computeProbability` switches to the generic algorithm while
PDF/CDF keep using Poisson's summation formula regardless of the size. Is
there any observable discontinuity at the threshold?

**Method.** :file:`p4_smallsize.py` sums N uniform U(0,1) atoms for
N = 99, 100, 101 and compares CDF values against a 10^6 sample Monte Carlo
reference, the local consistency :math:`p(x)\cdot 2h` versus
:math:`P([x-h, x+h])`, and the timings.

**Results.**

=====  ===================  ============  ============  ========
N      cdf error vs MC      local cons.   t(prob)       t(pdf)
=====  ===================  ============  ============  ========
99     -3.3e-05             2.74e-11      0.003 s       0.000 s
100    -3.1e-04             2.67e-11      0.001 s       0.001 s
101    +5.6e-04             2.61e-11      0.001 s       0.001 s
=====  ===================  ============  ============  ========

All deviations are within the Monte Carlo noise (standard deviation
~3.7e-4).

**Conclusion / recommendation.** No discontinuity, no cost jump: the
asymmetry is benign and mainly needs a documentation mention (now covered
by the ``SmallSize`` entry of the class Notes).


5. Remaining simplification opportunities
-----------------------------------------

5a. Opposite-weight Poisson atoms into a Skellam atom
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Currently ``Poi(l1) - Poi(l2)`` (unit weights) fuses into a truncated
FiniteDiscreteDistribution (43 points for l1=2, l2=3, tails beyond the
truncation dropped). Since
:math:`\sum_i Poi(\lambda_i) - \sum_j Poi(\mu_j) = Skellam(\sum\lambda,
\sum\mu)` **when all weights are +/-1**, substituting a single
:class:`~openturns.Skellam` atom would be exact and analytical.

:file:`p5_skellam.py` verifies the equivalence: mean/variance match
(-1/5) exactly and the PMFs agree to better than 1e-12 on [-3, 8).
Benefits: analytical status of the resulting distribution (fast paths),
no support truncation, fewer atoms. Restriction: invalid as soon as a
weight differs from +/-1 (scaling a Poisson leaves the family).

*Post-scriptum.* This substitution is now implemented: opposite-weight
Poisson groups fuse into a single Skellam atom weighted by the absolute
value of their weights.

5b. Cost of the discrete/continuous pairing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When both continuous and discrete atoms exist, each pair is replaced by a
Mixture of translated copies of the continuous atom (one component per
support point). :file:`p5_pairing.py` compares pointwise evaluations
against the unfused mixture (``SimplifyAtoms=False``):

====================================  =======  =============  ===========
configuration                        atoms    pdf ms/call    cdf ms/call
====================================  =======  =============  ===========
pairing ON (Binomial(2, 0.5))         2        1.1            0.1
pairing OFF                           3        < 0.1          < 0.1
pairing ON (Binomial(20, 0.5))        2        5.4            1.0
====================================  =======  =============  ===========

For this representative configuration the pairing is a *pessimization* for
pointwise evaluations: the generic Poisson summation on 3 raw atoms is below the timer
resolution and beats a
Mixture whose component count grows with the discrete support size.

**Recommendation.** Gate the pairing on a cost model (e.g. skip it when the
total number of mixture components exceeds a few units or when the
continuous atoms have cheap characteristic functions), or benchmark before
enabling; the current unconditional rule may slow down typical
mixed discrete/continuous models.


Reproducing
-----------

Each script next to this file is standalone:

.. code-block:: bash

    python3 p1_fusion.py     # fusion support sizes + MaximumSupportSize sweep
    python3 p2_tails.py      # discarded tail mass for heavy-tailed atoms
    python3 p3_cache.py      # timing of pointwise PDF vs DefaultMaxSize
    python3 p4_smallsize.py  # accuracy/cost across the SmallSize threshold
    python3 p5_skellam.py    # Skellam substitution equivalence check
    python3 p5_pairing.py    # cost of the discrete/continuous pairing

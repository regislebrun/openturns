%feature("docstring") OT::SubsampleChristoffelExperiment
R"RAW(Christoffel subsample experiment.

.. warning::
    This class is experimental and likely to be modified in future releases.
    To use it, import the ``openturns.experimental`` submodule.

Given an orthonormal basis of :math:`L^2(D, \mu)` and a target size
:math:`n`, the experiment draws a pool of points from the
:class:`~openturns.experimental.ChristoffelDistribution` associated to the
first :math:`m` basis functions, where :math:`m` is deduced from

.. math::
    n = C \, m \, \log m,

with :math:`C = 10` by default, then greedily thins the pool down to
:math:`n` points. The default ``Removal`` method drops pool points while
keeping the largest smallest Gramian eigenvalue and no reweighting; the
``Barrier`` method is a forward barrier greedy selection with reweighting,
clamped strictly inside the current spectrum so rank-deficient steps stay
well-defined. Each kept point :math:`x_i` carries a weight proportional
to :math:`m / k_m(x_i)`, the density ratio of the reference measure
over the Christoffel distribution, times the barrier weight for the
``Barrier`` method.

Available constructors:
    SubsampleChristoffelExperiment(*basis, size*)

Parameters
----------
basis : :class:`~openturns.OrthogonalBasis`
    Orthonormal basis defining the approximation space. Its embedded
    measure :math:`\mu` must be continuous.
size : positive int
    Number :math:`n` of points of the experiment.

See Also
--------
openturns.WeightedExperiment

Notes
-----
The following :class:`~openturns.ResourceMap` keys are used:

- ``SubsampleChristoffelExperiment-SamplingFactor`` (``Scalar``, default: ``10.0``): constant :math:`C` of the sizing rule.
- ``SubsampleChristoffelExperiment-PoolOversamplingFactor`` (``Scalar``, default: ``2.0``): candidate pool size as a multiple of the target size.
- ``SubsampleChristoffelExperiment-FrameTolerance`` (``Scalar``, default: ``0.5``): half-width of the accepted frame eigenvalue band around 1.
- ``SubsampleChristoffelExperiment-ThinningMethod`` (``String``, default: ``Removal``): greedy thinning method. Possible values: ``Barrier``, ``Removal``.
- ``SubsampleChristoffelExperiment-BarrierStep`` (``Scalar``, default: ``1.0``): barrier advance per selection step of the ``Barrier`` method.
- ``SubsampleChristoffelExperiment-BarrierRegularization`` (``Scalar``, default: ``1e-8``): safety gap keeping the barriers strictly inside the spectrum.

Examples
--------
>>> import openturns as ot
>>> import openturns.experimental as otexp
>>> ot.RandomGenerator.SetSeed(0)
>>> factory = ot.OrthogonalProductPolynomialFactory([ot.Uniform(-1.0, 1.0)])
>>> basis = ot.OrthogonalBasis(factory)
>>> experiment = otexp.SubsampleChristoffelExperiment(basis, 60)
>>> print(experiment.getSpaceDimension())
4
>>> sample, weights = experiment.generateWithWeights()
>>> print(len(sample), len(weights))
60 60)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::SubsampleChristoffelExperiment::getBasis
"Accessor to the orthogonal basis.

Returns
-------
basis : :class:`~openturns.OrthogonalBasis`
    Orthogonal basis defining the approximation space."

// ---------------------------------------------------------------------

%feature("docstring") OT::SubsampleChristoffelExperiment::getSpaceDimension
"Accessor to the space dimension deduced from the target size.

Returns
-------
spaceDimension : positive int
    Number of leading basis functions used."

// ---------------------------------------------------------------------

%feature("docstring") OT::SubsampleChristoffelExperiment::DeduceSpaceDimension
"Deduce the space dimension from a target size.

Parameters
----------
size : positive int
    Target number of points.

Returns
-------
spaceDimension : positive int
    Largest integer with factor times spaceDimension times
    log(spaceDimension) not larger than size."

// ---------------------------------------------------------------------

%feature("docstring") OT::SubsampleChristoffelExperiment::ComputeGramian
"Gramian of a subset of design points.

Parameters
----------
features : 2-d sequence of float
    Basis function values, one row per pool point.
weights : sequence of float
    Weight of each pool point.
kept : sequence of int
    Indices of the points kept in the subset.

Returns
-------
gramian : :class:`~openturns.SymmetricMatrix`
    Mean over the kept points of weight times outer product."

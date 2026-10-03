%feature("docstring") OT::ChristoffelDistribution
R"RAW(Christoffel distribution.

.. warning::
    This class is experimental and likely to be modified in future releases.
    To use it, import the ``openturns.experimental`` submodule.

The Christoffel distribution associated to the first :math:`m` functions
:math:`(L_1, \dots, L_m)` of an orthonormal basis of :math:`L^2(D, \mu)`
is the probability measure :math:`\sigma` defined by

.. math::
    \mathrm{d}\sigma(x) = \frac{k_m(x)}{m} \, \mathrm{d}\mu(x),
    \quad k_m(x) = \sum_{j=1}^m L_j(x)^2,

where :math:`k_m` is the Christoffel function and :math:`\mu` is the
reference measure embedded in the basis. When :math:`\mu` is absolutely
continuous with density :math:`p_\mu`, the density of :math:`\sigma` is

.. math::
    p(x) = p_\mu(x) \, \frac{k_m(x)}{m}.

Available constructors:
    ChristoffelDistribution(*basis, size*)

Parameters
----------
basis : :class:`~openturns.OrthogonalBasis`
    Orthonormal basis defining the approximation space. Its embedded
    measure :math:`\mu` must be continuous.
size : positive int
    Dimension :math:`m` of the approximation space, ie the number of
    leading basis functions used.

See Also
--------
openturns.Distribution

Notes
-----
The following :class:`~openturns.ResourceMap` keys are used:

- ``ChristoffelDistribution-OptimizationAlgorithm`` (``String``, default: ``Cobyla``): optimization algorithm used by the internal ratio-of-uniforms sampler.
- ``ChristoffelDistribution-RatioUniformCandidateNumber`` (``UnsignedInteger``, default: ``10000``): number of candidate points of the internal ratio-of-uniforms sampler.
- ``ChristoffelDistribution-RatioUniformMaxDimension`` (``UnsignedInteger``, default: ``5``): input dimension above which sampling falls back to crude rejection instead of ratio-of-uniforms.
- ``ChristoffelDistribution-KnSamplingSize`` (``UnsignedInteger``, default: ``100000``): sample size of the Monte-Carlo estimate of the stability factor.
- ``ChristoffelDistribution-KnSafetyFactor`` (``Scalar``, default: ``1.0``): multiplicative safety margin of the stability factor estimate.
- ``ChristoffelDistribution-SliceGridSize`` (``UnsignedInteger``, default: ``1000``): grid size of the univariate slice maxima in sequential conditional sampling.

Examples
--------
>>> import openturns as ot
>>> import openturns.experimental as otexp
>>> factory = ot.OrthogonalProductPolynomialFactory([ot.Uniform(-1.0, 1.0)])
>>> basis = ot.OrthogonalBasis(factory)
>>> distribution = otexp.ChristoffelDistribution(basis, 3)
>>> print(distribution.getSize())
3
>>> print(f"{distribution.computePDF([0.5]):.6f}")
0.304688)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::getOrthogonalBasis
"Accessor to the orthogonal basis.

Returns
-------
basis : :class:`~openturns.OrthogonalBasis`
    Orthogonal basis defining the approximation space."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::getSize
"Accessor to the space dimension.

Returns
-------
size : positive int
    Number of leading basis functions used."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::getMeasure
"Accessor to the reference measure embedded in the basis.

Returns
-------
measure : :class:`~openturns.Distribution`
    Reference measure with respect to which the basis is orthonormal."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::getBasis
"Accessor to the finite basis of leading functions.

Returns
-------
basis : :class:`~openturns.Basis`
    First functions of the orthogonal basis."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::computeChristoffel
"Evaluate the Christoffel function.

Parameters
----------
point : sequence of float or 2-d sequence of float
    Point or sample where the Christoffel function is evaluated.

Returns
-------
values : float or :class:`~openturns.Sample`
    Sum of the squared basis functions at the given point or sample."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::computeKn
"Monte-Carlo estimate of the stability factor.

Returns
-------
kn : float
    Maximum of the Christoffel function over a sample drawn from
    the reference measure, times the safety margin."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::computeConditionalPDF
"Conditional PDF of a component given the previous ones.

Parameters
----------
x : float
    Value of the conditioned component.
y : sequence of float
    Values of the conditioning components.

Returns
-------
value : float
    Closed slice form in the tensor-product case, generic
    implementation otherwise."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::computeConditionalCDF
"Conditional CDF of a component given the previous ones.

Parameters
----------
x : float
    Value of the conditioned component.
y : sequence of float
    Values of the conditioning components.

Returns
-------
value : float
    Slice quadrature in the tensor-product case, generic
    implementation otherwise."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::computeSequentialConditionalPDF
"Conditional PDFs of each component given the previous ones.

Parameters
----------
x : sequence of float
    Point where the sequential conditionals are evaluated.

Returns
-------
values : :class:`~openturns.Point`
    Fan-out to the scalar conditional PDFs in the tensor-product
    case, generic implementation otherwise."

// ---------------------------------------------------------------------

%feature("docstring") OT::ChristoffelDistribution::computeSequentialConditionalCDF
"Conditional CDFs of each component given the previous ones.

Parameters
----------
x : sequence of float
    Point where the sequential conditionals are evaluated.

Returns
-------
values : :class:`~openturns.Point`
    Fan-out to the scalar conditional CDFs in the tensor-product
    case, generic implementation otherwise."

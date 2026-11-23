%feature("docstring") OT::TruncatedDistribution
R"RAW(Truncated distribution.

Assuming a random variable *X* is truncated with two bounds *a* and *b*,
its probability density function is defined as:

.. math::

      \forall y \in \Rset, p_Y(y) =
      \begin{array}{|ll}
        0 &  \mbox{for } y \geq b  \mbox{ or }  y \leq a\\
        \displaystyle \frac{1}{F_X(b) - F_X(a)}\, p_X(y) & \mbox{for } y\in[a,b]
      \end{array}

and its cumulative distribution function is defined as:

.. math::

    \forall y \in \Rset, F_Y(y) = \Prob{X<y\, / \, X\in[a,b]} =
    \begin{array}{|ll}
      1 & \mbox{for } y \geq b, \\
      0 & \mbox{for } y \leq a, \\
      \displaystyle \frac{F_X(y) - F_X(a)}{F_X(b) - F_X(a)} & \mbox{for } y\in[a,b]
    \end{array}

Available constructors:
    TruncatedDistribution(*distribution, bound, side, threshold*)

    TruncatedDistribution(*distribution, truncationInterval, threshold*)

Parameters
----------
distribution : :class:`~openturns.Distribution`
    The distribution to be truncated.
bound : float
    Scalar used when one of the bounds of the distribution is not finite.

    If *side = ot.TruncatedDistribution.LOWER*: the interval
    :math:`[bound, +\infty)` is used to bound the distribution.

    If *side = ot.TruncatedDistribution.UPPER*: the interval
    :math:`(-\infty, bound]` is used to bound the distribution.
side : int, optional
    Integer used when one of the bounds of the distribution is not finite.
    (Amongst ot.TruncatedDistribution.LOWER and ot.TruncatedDistribution.UPPER).
    By default, equal to ot.TruncatedDistribution.LOWER.
truncationInterval : :class:`~openturns.Interval`
    Define a new distribution range :math:`[lowerBound, upperBound]` with
    finite bounds.
threshold : float, :math:`\tau \in [0, 1]`
    This threshold is used when generating an evaluation of the distribution.
    By default, equal to 0.5, which is the value of the
    ``TruncatedDistribution-DefaultThresholdRealization`` :class:`~openturns.ResourceMap` key.
    We note :math:`F` the CDF of the non truncated distribution:

    If :math:`F(upperBound)-F(lowerBound)<\tau`: a CDF inversion is performed.

    If :math:`F(upperBound)-F(lowerBound)>\tau`: rejection.

Notes
-----
When the distribution to be truncated is a multivariate normal distribution
with correlated components, the mean and the covariance of the truncated
distribution have no closed form: they are estimated by the weighted Genz
quasi-Monte Carlo algorithm, in a single pass shared by both. The estimate is
exact in dimension 1, where the truncation simplifies to a
:class:`~openturns.TruncatedNormal`, and when the components are independent,
where it is the product of univariate truncated normal distributions.

The following :class:`~openturns.ResourceMap` keys are used:

- ``TruncatedDistribution-DefaultThresholdRealization`` (``Scalar``, default: ``0.5``): default value of the threshold :math:`\tau` used to choose between the rejection method and the CDF inversion.
- ``TruncatedDistribution-MinimumEffectiveSampleSize`` (``Scalar``, default: ``100``): the Genz estimate of the moments warns when the effective sample size :math:`\left(\sum_i p_i\right)^2 / \sum_i p_i^2` of the box weights :math:`p_i` falls below this value, as the relative error of the estimate is then of the order of the inverse square root of the effective sample size.

Examples
--------
>>> import openturns as ot
>>> distribution = ot.TruncatedDistribution(ot.Normal(2.0, 1.5), ot.Interval(1.0, 4.0))

Draw a sample:

>>> sample = distribution.getSample(5)
)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::getDistribution
"Get the distribution.

Returns
-------
distribution : :class:`~openturns.Distribution`
    The distribution to be truncated."

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::setDistribution
"Set the distribution.

Parameters
----------
distribution : :class:`~openturns.Distribution`
    The distribution to be truncated."

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::setBounds
"Set the bounds.

Parameters
----------
bounds : :class:`~openturns.Interval`
    Truncation bounds."

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::getBounds
"Set the bounds.

Returns
-------
bounds : :class:`~openturns.Interval`
    Truncation bounds."

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::getThresholdRealization
R"RAW(Get the threshold used when generating an evaluation of the distribution.

Returns
-------
threshold : float, :math:`\tau \in [0, 1]`
    This threshold is used when generating an evaluation of the distribution.
    We note :math:`F` the CDF of the non truncated distribution:

    If :math:`F(upperBound)-F(lowerBound)<\tau`: a CDF inversion is performed.

    If :math:`F(upperBound)-F(lowerBound)>\tau`: rejection.)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::setThresholdRealization
R"RAW(Set the threshold used when generating an evaluation of the distribution.

Parameters
----------
threshold : float, :math:`\tau \in [0, 1]`
    This threshold is used when generating an evaluation of the distribution.
    We note :math:`F` the CDF of the non truncated distribution:

    If :math:`F(upperBound)-F(lowerBound)<\tau`: a CDF inversion is performed.

    If :math:`F(upperBound)-F(lowerBound)>\tau`: rejection.)RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::TruncatedDistribution::getSimplifiedVersion
"Accessor to the simplified distribution.

Drops unneeded truncation if the distribution given as argument is already
truncated by nature (Uniform, etc).

Returns
-------
simplified : :class:`~openturns.Distribution`
    The simplified distribution."

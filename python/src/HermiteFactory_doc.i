%feature("docstring") OT::HermiteFactory
R"RAW(Hermite specific orthonormal univariate polynomial family.

For the :class:`~openturns.Normal` distribution.

Any sequence of orthogonal polynomials has a recurrence formula relating any
three consecutive polynomials as follows:

.. math::

    P_{i + 1} = (a_i x + b_i) P_i + c_i P_{i - 1}, \quad 1 < i

The recurrence coefficients for the Hermite polynomials come analytically and
read:

.. math::

    \begin{array}{rcl}
        a_i & = & \displaystyle \frac{1}{\sqrt{i + 1}} \\
        b_i & = & 0 \\
        c_i & = & \displaystyle - \sqrt{\frac{i}{i + 1}}
    \end{array}, \quad 1 < i

The nodes and weights of the associated Gauss-Hermite quadrature rule are
computed by the fast Hermite rule mapped to the measure: polished
eigensolver below 256 nodes, Townsend-Trogdon-Olver Airy expansion above.

See also
--------
UniVariateDistributionPolynomialFactory

Notes
-----
Above 256 nodes the rule reaches a relative accuracy better than
``5e-13`` and is faster than the generic solver; the switch threshold
comes from the ``fast_gauss`` benchmark. The following
:class:`~openturns.ResourceMap` key is used:

- ``FastHermite-AsymptoticThreshold`` (``UnsignedInteger``, default:
  ``256``): number of nodes from which the asymptotic expansion is used.
  Set it to a large number to force the generic polished eigensolver:
  improved accuracy beyond ``5e-13`` at the price of a much larger CPU
  effort.

Examples
--------
>>> import openturns as ot
>>> polynomial_factory = ot.HermiteFactory()
>>> for i in range(3):
...     print(polynomial_factory.build(i))
1
X
-0.707107 + 0.707107 * X^2

>>> polynomial_factory = ot.HermiteFactory(1.0, 2.0)
>>> print(polynomial_factory)
class=HermiteFactory measure=class=Normal name=Normal dimension=1 mean=class=Point name=Unnamed dimension=1 values=[1] sigma=class=Point name=Unnamed dimension=1 values=[2] correlationMatrix=class=CorrelationMatrix dimension=1 implementation=class=MatrixImplementation name=Unnamed rows=1 columns=1 values=[1])RAW"

%feature("docstring") OT::AbdoRackwitz
"Abdo-Rackwitz solver.

This solver uses first derivative information and can only be used to solve level function problems.

Available constructors:
    AbdoRackwitz(*problem*)

    AbdoRackwitz(*problem, tau, omega, smooth*)

Parameters
----------
problem : :class:`~openturns.OptimizationProblem`
    Optimization problem to solve.
tau : float
    Multiplicative decrease of linear step.
omega : float
    Armijo factor.
smooth : float
    Growing factor in penalization term.

Notes
-----
The algorithm stops when one of the two following pairs of criteria is
met, in which case it has converged:

- the absolute error and the relative error are both below their maximum,
  set through :meth:`setMaximumAbsoluteError` and
  :meth:`setMaximumRelativeError`;
- the residual error and the constraint error are both below their
  maximum, set through :meth:`setMaximumResidualError` and
  :meth:`setMaximumConstraintError`.

Because the rule is a disjunction, setting some of these thresholds to
zero only disables the corresponding pair: the other pair still applies
with its own default thresholds.

See also
--------
Cobyla, SQP, TNC, NLopt

Examples
--------
>>> import openturns as ot
>>> model = ot.SymbolicFunction(['E', 'F', 'L', 'I'], ['-F*L^3/(3*E*I)'])
>>> problem = ot.NearestPointProblem(model, 5.0)
>>> algo = ot.AbdoRackwitz(problem)
>>> algo.setStartingPoint([1.0] * 4)
>>> algo.run()
>>> result = algo.getResult()"

// ---------------------------------------------------------------------

%feature("docstring") OT::AbdoRackwitz::getTau
"Accessor to tau parameter.

Returns
-------
tau : float
    Multiplicative decrease of linear step."

// ---------------------------------------------------------------------

%feature("docstring") OT::AbdoRackwitz::setTau
"Accessor to tau parameter.

Parameters
----------
tau : float
    Multiplicative decrease of linear step."

// ---------------------------------------------------------------------

%feature("docstring") OT::AbdoRackwitz::getOmega
"Accessor to omega parameter.

Returns
-------
omega : float
    Armijo factor."

// ---------------------------------------------------------------------

%feature("docstring") OT::AbdoRackwitz::setOmega
"Accessor to omega parameter.

Parameters
----------
omega : float
    Armijo factor."

// ---------------------------------------------------------------------

%feature("docstring") OT::AbdoRackwitz::getSmooth
"Accessor to smooth parameter.

Returns
-------
smooth : float
    Growing factor in penalization term."

// ---------------------------------------------------------------------

%feature("docstring") OT::AbdoRackwitz::setSmooth
"Accessor to smooth parameter.

Parameters
----------
smooth : float
    Growing factor in penalization term."


%feature("docstring") OT::ComposedFunction
R"RAW(Composed function.

Available constructor:
    ComposedFunction(*f, g*)

The function is the composed function :math:`f\circ g`.

Parameters
----------
f,g : two :class:`~openturns.Function`
    The functions to compose.

Examples
--------
>>> import openturns as ot
>>> g = ot.SymbolicFunction(['x1', 'x2'],
...                           ['x1 + x2','3 * x1 * x2'])
>>> f = ot.SymbolicFunction(['x1', 'x2'], ['2 * x1 - x2'])
>>> composed = ot.ComposedFunction(f, g)
>>> print(composed([3, 4]))
[-22])RAW"

// ---------------------------------------------------------------------

%feature("docstring") OT::ComposedFunction::getLeftFunction
"Left function accessor.

Returns
-------
leftFunction : :class:`~openturns.Function`
    The outer function :math:`f` of the composition :math:`f \\circ g`."

// ---------------------------------------------------------------------

%feature("docstring") OT::ComposedFunction::setLeftFunction
"Set the left function, ie the outer function of the composition.

Parameters
----------
leftFunction : :class:`~openturns.Function`
    The outer function of the composition."

// ---------------------------------------------------------------------

%feature("docstring") OT::ComposedFunction::getRightFunction
"Right function accessor.

Returns
-------
rightFunction : :class:`~openturns.Function`
    The inner function :math:`g` of the composition :math:`f \\circ g`."

// ---------------------------------------------------------------------

%feature("docstring") OT::ComposedFunction::setRightFunction
"Set the right function, ie the inner function of the composition.

Parameters
----------
rightFunction : :class:`~openturns.Function`
    The inner function of the composition."

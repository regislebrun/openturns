#! /usr/bin/env python

import openturns as ot
import openturns.testing as ott
import math as m

vertices = ot.Sample([[0.0], [0.1], [0.2], [0.3]])
size = vertices.getSize()


class TestAssembly:
    def __init__(self, vertices):
        self.vertices = vertices

    def __call__(self, i, j):
        diff = self.vertices[i][0] - self.vertices[j][0]
        return m.exp(-abs(diff) / 0.1)


factory = ot.HODLRMatrixFactory()
parameters = ot.HODLRMatrixParameters()
parameters.setAssemblyEpsilon(1.0e-6)
parameters.setMinLeafSize(2)

hodlr = factory.build(vertices, 1, True, parameters)
print("built n=", hodlr.getNbRows())

assemblyFunc = TestAssembly(vertices)
hodlr.assembleReal(assemblyFunc, 'L')
print("assembled")

print("about to factorize...")
hodlr.factorize('LU')
print("factorized, logDet=", hodlr.logDeterminant())

b = ot.Point(size, 1.0)
hodlr2 = factory.build(vertices, 1, True, parameters)
hodlr2.assembleReal(assemblyFunc, 'L')
with ott.assert_raises(TypeError):
    hodlr2.solve(b)
